/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package speclib.io;

import java.io.IOException;
import java.sql.*;
import java.util.*;

public class ParquetToSpecLib {
    
    private final String parquetFilePath;
    private final Map<String, String> isoformToGeneMap;
    private final int version;
    private final boolean genDecoys;
    private final boolean genCharges;
    private final boolean inferProteotypicity;
    
    private static class TempPrecursorData {
        List<Product> fragments;
        float precursorMz;
        String proteinGroup;
        String peptideSequence;
        String modifiedPeptideSequence;
        short precursorCharge;
        float normalizedRetentionTime;
        float precursorIonMobility;
        short proteotypic;
    }

    public ParquetToSpecLib(String parquetFilePath, Map<String, String> isoformToGeneMap, int version, boolean genDecoys, boolean genCharges, boolean inferProteotypicity) {
        this.parquetFilePath = parquetFilePath;
        this.isoformToGeneMap = isoformToGeneMap != null ? isoformToGeneMap : Collections.emptyMap();
        this.version = version;
        this.genDecoys = genDecoys;
        this.genCharges = genCharges;
        this.inferProteotypicity = inferProteotypicity;
    }

    public int getVersion() {
        return version;
    }
    
    private static boolean isSwissprot(String isoformId) {
        if (isoformId == null || isoformId.isEmpty()) {
            return false;
        }
        return isoformId.startsWith("sp|");
    }
    
    public SpectralLibrary convert() throws IOException {
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:")) {            
            long rowCount = getRowCount(conn);
            return buildSpectralLibraryStreaming(conn, rowCount);
        } catch (SQLException e) {
            throw new IOException("Error reading Parquet file with DuckDB: " + e.getMessage(), e);
        }
    }
    
    private long getRowCount(Connection conn) throws SQLException {
        String countQuery = "SELECT COUNT(*) FROM read_parquet(?)";
        try (PreparedStatement stmt = conn.prepareStatement(countQuery)) {
            stmt.setString(1, parquetFilePath);
            try (ResultSet rs = stmt.executeQuery()) {
                if (rs.next()) {
                    return rs.getLong(1);
                }
            }
        }
        return 1L << 15;
    }

    private SpectralLibrary buildSpectralLibraryStreaming(Connection conn, long rowCount) throws SQLException, IOException {
        SpectralLibrary library = new SpectralLibrary();
        
        int estimatedPrecursors = (int) (rowCount / 16);
        
        Set<String> uniqueIsoforms = new TreeSet<>();
        Set<String> uniqueProteinGroups = new TreeSet<>();
        Set<String> uniqueGenes = new TreeSet<>();
        
        List<String> precursors = new ArrayList<>(estimatedPrecursors);
        List<LibraryEntry> entries = new ArrayList<>(estimatedPrecursors);
        List<TempPrecursorData> tempData = new ArrayList<>(estimatedPrecursors);
        
        double minRT = Double.MAX_VALUE;
        double maxRT = Double.MIN_VALUE;
        
        String query = "SELECT " +
            "ModifiedPeptideSequence, " +
            "PrecursorCharge, " +
            "FIRST(PrecursorMz), " +
            "FIRST(ProteinId), " +
            "FIRST(PeptideSequence), " +
            "FIRST(NormalizedRetentionTime), " +
            "FIRST(COALESCE(TRY_CAST(PrecursorIonMobility AS FLOAT), 0.0)), " +
            "FIRST(Proteotypic), " +
            "list(ProductMz), " +
            "list(LibraryIntensity), " +
            "list(FragmentCharge), " +
            "list(CASE " +
            "  WHEN lower(FragmentType) = 'b' THEN 1 " +
            "  WHEN lower(FragmentType) = 'y' THEN 2 " +
            "  ELSE 0 " +
            "END), " +
            "list(FragmentSeriesNumber), " +
            "list(CASE " +
            "  WHEN FragmentLossType IS NULL OR FragmentLossType = '' OR lower(FragmentLossType) = 'noloss' THEN 0 " +
            "  WHEN lower(FragmentLossType) IN ('h2o', '-h2o') THEN 1 " +
            "  WHEN lower(FragmentLossType) IN ('nh3', '-nh3') THEN 2 " +
            "  WHEN lower(FragmentLossType) IN ('co', '-co') THEN 3 " +
            "  WHEN lower(FragmentLossType) IN ('n', '-n') THEN 4 " +
            "  WHEN lower(FragmentLossType) IN ('other', '-other') THEN 5 " +
            "  ELSE 0 " +
            "END) " +
            "FROM read_parquet(?) " +
            "GROUP BY ModifiedPeptideSequence, PrecursorCharge " +
            "ORDER BY ModifiedPeptideSequence, PrecursorCharge";
        
        try (PreparedStatement stmt = conn.prepareStatement(query,
                ResultSet.TYPE_FORWARD_ONLY, 
                ResultSet.CONCUR_READ_ONLY)) {
            stmt.setFetchSize(1 << 10);
            stmt.setString(1, parquetFilePath);
            
            try (ResultSet rs = stmt.executeQuery()) {
                while (rs.next()) {
                    String modifiedPeptideSequence = rs.getString(1);
                    short precursorCharge = rs.getShort(2);
                    float precursorMz = rs.getFloat(3);
                    String proteinGroup = rs.getString(4);
                    String peptideSequence = rs.getString(5);
                    float normalizedRetentionTime = rs.getFloat(6);
                    float precursorIonMobility = rs.getFloat(7);
                    short proteotypic = rs.getShort(8);
                    
                    Array productMzArray = rs.getArray(9);
                    Array libraryIntensityArray = rs.getArray(10);
                    Array fragmentChargeArray = rs.getArray(11);
                    Array fragmentTypeArray = rs.getArray(12);
                    Array fragmentSeriesNumberArray = rs.getArray(13);
                    Array fragmentLossTypeArray = rs.getArray(14);
                    
                    Object[] mzs = (Object[]) productMzArray.getArray();
                    Object[] intensities = (Object[]) libraryIntensityArray.getArray();
                    Object[] charges = (Object[]) fragmentChargeArray.getArray();
                    Object[] types = (Object[]) fragmentTypeArray.getArray();
                    Object[] series = (Object[]) fragmentSeriesNumberArray.getArray();
                    Object[] losses = (Object[]) fragmentLossTypeArray.getArray();
                    
                    int fragmentCount = mzs.length;
                    List<Product> fragments = new ArrayList<>(fragmentCount);
                    for (int i = 0; i < fragmentCount; i++) {
                        byte type = ((Number)types[i]).byteValue();
                        if (type == 0) throw new IOException("Unsupported fragment type found in parquet file.");
                        
                        fragments.add(new Product(
                            ((Number)mzs[i]).floatValue(),
                            ((Number)intensities[i]).floatValue(),
                            ((Number)charges[i]).byteValue(),
                            type,
                            ((Number)series[i]).byteValue(),
                            ((Number)losses[i]).byteValue()
                        ));
                    }
                    
                    if (!fragments.isEmpty()) {
                        TempPrecursorData data = new TempPrecursorData();
                        data.fragments = fragments;
                        data.precursorMz = precursorMz;
                        data.proteinGroup = proteinGroup != null ? proteinGroup : "";
                        data.peptideSequence = peptideSequence;
                        data.modifiedPeptideSequence = modifiedPeptideSequence;
                        data.precursorCharge = precursorCharge;
                        data.normalizedRetentionTime = normalizedRetentionTime;
                        data.precursorIonMobility = precursorIonMobility;
                        data.proteotypic = proteotypic;
                        
                        tempData.add(data);

                        for (String isoform : data.proteinGroup.split(";")) {
                            uniqueIsoforms.add(isoform);
                            uniqueGenes.add(isoformToGeneMap.getOrDefault(isoform, isoform));
                        }
                        uniqueProteinGroups.add(data.proteinGroup);
                    }
                    
                    if (normalizedRetentionTime < minRT) minRT = normalizedRetentionTime;
                    if (normalizedRetentionTime > maxRT) maxRT = normalizedRetentionTime;
                }
            }
        }
        
        List<String> sortedIsoforms = new ArrayList<>(uniqueIsoforms);
        Map<String, Integer> isoformIdxMap = new HashMap<>(sortedIsoforms.size());
        for (int i = 0; i < sortedIsoforms.size(); ++i) {
            isoformIdxMap.put(sortedIsoforms.get(i), i);
        }

        List<String> sortedProteinGroups = new ArrayList<>(uniqueProteinGroups);
        Map<String, Integer> proteinGroupIdxMap = new HashMap<>(sortedProteinGroups.size());
        for (int i = 0; i < sortedProteinGroups.size(); ++i) {
            proteinGroupIdxMap.put(sortedProteinGroups.get(i), i);
        }
        
        List<String> sortedGenes = new ArrayList<>(uniqueGenes);
        Map<String, Integer> geneIdxMap = new HashMap<>(sortedGenes.size());
        for (int i = 0; i < sortedGenes.size(); ++i) {
            geneIdxMap.put(sortedGenes.get(i), i);
        }

        Map<String, Set<Integer>> isoformToPrecursors = new HashMap<>();
        for (int precursorIndex = 0; precursorIndex < tempData.size(); ++precursorIndex) {
            TempPrecursorData data = tempData.get(precursorIndex);
            for (String isoform : data.proteinGroup.split(";")) {
                isoformToPrecursors.computeIfAbsent(isoform, k -> new LinkedHashSet<>()).add(precursorIndex);
            }
        }

        List<Isoform> isoforms = new ArrayList<>(sortedIsoforms.size());
        for (String s : sortedIsoforms) {
            String gene = isoformToGeneMap.getOrDefault(s, s);
            Isoform isoform = new Isoform();
            isoform.setId(s);
            isoform.setName("");
            isoform.setGene(gene);
            isoform.setDescription("");
            isoform.setNameIndex(0);
            isoform.setGeneIndex(geneIdxMap.get(gene));
            isoform.setSwissprot(isSwissprot(s));
            isoform.setPrecursors(isoformToPrecursors.getOrDefault(s, new LinkedHashSet<>()));
            isoforms.add(isoform);
        }
        
        Map<Integer, List<Integer>> proteinGroupToPrecursors = new HashMap<>(sortedProteinGroups.size());
        for (int precursorIndex = 0; precursorIndex < tempData.size(); ++precursorIndex) {
            TempPrecursorData data = tempData.get(precursorIndex);
            int proteinGroupIdx = proteinGroupIdxMap.get(data.proteinGroup);
            
            proteinGroupToPrecursors.computeIfAbsent(proteinGroupIdx, k -> new ArrayList<>()).add(precursorIndex);
            
            String precursorId = data.modifiedPeptideSequence + data.precursorCharge;
            precursors.add(precursorId);
            
            LibraryEntry libEntry = new LibraryEntry();
            libEntry.setName(precursorId);
            libEntry.setPgIndex(proteinGroupIdx);
            libEntry.setProteotypic(data.proteotypic);
            
            Precursor precursor = new Precursor();
            precursor.setIndex(precursorIndex);
            precursor.setCharge(data.precursorCharge);
            precursor.setLength(data.peptideSequence != null ? data.peptideSequence.length() : 0);
            precursor.setMz(data.precursorMz);
            precursor.setiRT(data.normalizedRetentionTime);
            precursor.setiIM(data.precursorIonMobility);
            precursor.setFragments(data.fragments);
            
            libEntry.setTarget(precursor);
            entries.add(libEntry);
        }

        List<ProteinGroup> proteinGroups = new ArrayList<>(sortedProteinGroups.size());
        for (int i = 0; i < sortedProteinGroups.size(); ++i) {
            String proteinGroupStr = sortedProteinGroups.get(i);
            List<Integer> precursorList = proteinGroupToPrecursors.get(i);

            List<Integer> isoformIndices = new ArrayList<>(1);
            List<String> genes = new ArrayList<>(1);
            List<Integer> geneIndices = new ArrayList<>(1);
            for (String isoform : proteinGroupStr.split(";")) {
                isoformIndices.add(isoformIdxMap.get(isoform));
                String gene = isoformToGeneMap.getOrDefault(isoform, isoform);
                genes.add(gene);
                geneIndices.add(geneIdxMap.get(gene));
            }
            
            ProteinGroup proteinGroup = new ProteinGroup();
            proteinGroup.setIds(proteinGroupStr);
            proteinGroup.setNames("");
            proteinGroup.setGenes(String.join(";", genes));
            proteinGroup.setPrecursors(precursorList != null ? precursorList : new ArrayList<>());
            proteinGroup.setIsoforms(isoformIndices);
            proteinGroup.setNameIndices(new ArrayList<>(0));
            proteinGroup.setGeneIndices(geneIndices);
            
            proteinGroups.add(proteinGroup);
        }

        library.setName("Converted from Parquet");
        library.setFastaNames("");
        library.setIsoforms(isoforms);
        library.setProteinGroups(proteinGroups);
        library.setPrecursors(precursors);
        library.setNames(new ArrayList<>(0));
        library.setGenes(sortedGenes);
        library.setiRTMin(minRT);
        library.setiRTMax(maxRT);
        library.setEntries(entries);
        library.setGenDecoys(genDecoys);
        library.setGenCharges(genCharges);
        library.setInferProteotypicity(inferProteotypicity);
        
        return library;
    }

    public void convertAndWrite(String outputPath) throws IOException {
        SpectralLibrary library = convert();
        
        long startWrite = System.nanoTime();
        DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
        writer.write(outputPath, version);
        long endWrite = System.nanoTime();
        System.out.println("Write time: " + (endWrite - startWrite) / 1_000_000 + " ms");
    }

    public static void main(String[] args) {
        String parquetFilePath = "I:\\test_transfer_learning\\dia_1\\fragpipe-predicted-speclib.parquet";
        String outputFilePath = "I:\\test_transfer_learning\\dia_1\\fragpipe-predicted-speclib.speclib";

        long startTime = System.nanoTime();

        ParquetToSpecLib converter = new ParquetToSpecLib(parquetFilePath, Collections.emptyMap(), -1, true, true, true);
        try {
            long convertStart = System.nanoTime();
            SpectralLibrary library = converter.convert();
            long convertEnd = System.nanoTime();
            System.out.println("Conversion time: " + (convertEnd - convertStart) / 1_000_000 + " ms");
            
            long writeStart = System.nanoTime();
            DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
            writer.write(outputFilePath, converter.getVersion());
            long writeEnd = System.nanoTime();
            System.out.println("Final write time: " + (writeEnd - writeStart) / 1_000_000 + " ms");
            
            System.out.println("Library stats: " + library.getEntries().size() + " precursors, " + 
                             library.getIsoforms().size() + " proteins");
        } catch (IOException e) {
            System.err.println("Error converting Parquet file: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
        long endTime = System.nanoTime();
        System.out.println("Conversion complete. Output written to " + outputFilePath);
        System.out.println("Total execution time: " + (endTime - startTime) / 1_000_000 + " ms");
    }
}
