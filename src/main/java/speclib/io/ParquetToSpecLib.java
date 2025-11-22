/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */
package speclib.io;

import java.io.IOException;
import java.sql.*;
import java.util.*;

public class ParquetToSpecLib {
    
    private final String parquetFilePath;
    private final Map<String, String> proteinToGeneMap;
    
    public ParquetToSpecLib(String parquetFilePath) {
        this(parquetFilePath, Collections.emptyMap());
    }

    public ParquetToSpecLib(String parquetFilePath, Map<String, String> proteinToGeneMap) {
        this.parquetFilePath = parquetFilePath;
        this.proteinToGeneMap = proteinToGeneMap != null ? proteinToGeneMap : Collections.emptyMap();
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
        int estimatedProteins = Math.max(1 << 13, (int) (rowCount / 128));
        
        Map<String, Integer> proteinIdMap = new LinkedHashMap<>(estimatedProteins);
        Map<String, Integer> geneNameMap = new LinkedHashMap<>(estimatedProteins);
        Map<Integer, List<Integer>> proteinToPrecursors = new LinkedHashMap<>(estimatedProteins);
        
        List<Isoform> proteins = new ArrayList<>(estimatedProteins);
        List<ProteinGroup> proteinIds = new ArrayList<>(estimatedProteins);
        List<String> precursors = new ArrayList<>(estimatedPrecursors);
        List<String> names = new ArrayList<>(estimatedProteins);
        List<String> genes = new ArrayList<>(estimatedProteins);
        List<LibraryEntry> entries = new ArrayList<>(estimatedPrecursors);
        
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
                int precursorIndex = 0;
                
                while (rs.next()) {
                    String modifiedPeptideSequence = rs.getString(1);
                    short precursorCharge = rs.getShort(2);
                    float precursorMz = rs.getFloat(3);
                    String proteinId = rs.getString(4);
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
                    
                    processPrecursor(
                        fragments,
                        precursorMz, proteinId, peptideSequence, modifiedPeptideSequence,
                        precursorCharge, normalizedRetentionTime, precursorIonMobility, proteotypic,
                        proteinIdMap, geneNameMap, proteinToPrecursors, proteins, names, genes, precursors, entries,
                        precursorIndex
                    );
                    
                    if (normalizedRetentionTime < minRT) minRT = normalizedRetentionTime;
                    if (normalizedRetentionTime > maxRT) maxRT = normalizedRetentionTime;
                    
                    precursorIndex++;
                }
            }
        }
        for (int i = 0; i < proteins.size(); i++) {
            Isoform protein = proteins.get(i);
            String geneName = protein.getGene();
            
            protein.setNameIndex(i < names.size() ? i : 0);
            protein.setGeneIndex(geneNameMap.getOrDefault(geneName, 0));
            
            List<Integer> precursorList = proteinToPrecursors.get(i);
            if (precursorList != null) {
                protein.setPrecursors(new LinkedHashSet<>(precursorList));
            }
            
            ProteinGroup pg = new ProteinGroup();
            pg.setIds(protein.getId());
            pg.setNames(protein.getName());
            pg.setGenes(protein.getGene());
            pg.setPrecursors(precursorList != null ? precursorList : new ArrayList<>());
            
            Set<Integer> proteinSet = new LinkedHashSet<>();
            proteinSet.add(i);
            pg.setProteins(proteinSet);
            
            List<Integer> nameIndices = new ArrayList<>(1);
            nameIndices.add(protein.getNameIndex());
            pg.setNameIndices(nameIndices);
            
            List<Integer> geneIndices = new ArrayList<>(1);
            geneIndices.add(protein.getGeneIndex());
            pg.setGeneIndices(geneIndices);
            
            proteinIds.add(pg);
        }
        
        if (names.isEmpty()) {
            for (Isoform protein : proteins) {
                names.add(protein.getName());
            }
        }

        library.setName("Converted from Parquet");
        library.setFastaNames("");
        library.setProteins(proteins);
        library.setProteinIds(proteinIds);
        library.setPrecursors(precursors);
        library.setNames(names);
        library.setGenes(genes);
        library.setiRTMin(minRT);
        library.setiRTMax(maxRT);
        library.setEntries(entries);
        library.setGenDecoys(false);
        library.setGenCharges(false);
        library.setInferProteotypicity(false);
        
        return library;
    }

    private void processPrecursor(
            List<Product> fragments,
            float precursorMz, String proteinId, String peptideSequence, String modifiedPeptideSequence,
            short precursorCharge, float normalizedRetentionTime, float precursorIonMobility, short proteotypic,
            Map<String, Integer> proteinIdMap, Map<String, Integer> geneNameMap, Map<Integer, List<Integer>> proteinToPrecursors,
            List<Isoform> proteins, List<String> names, List<String> genes, List<String> precursors, List<LibraryEntry> entries,
            int precursorIndex) {
        
        if (fragments.isEmpty()) {
            return;
        }
        
        String safeProteinId = proteinId != null ? proteinId : "";
        String geneName = proteinToGeneMap.getOrDefault(safeProteinId, safeProteinId);
        
        int proteinIdx;
        if (!proteinIdMap.containsKey(safeProteinId)) {
            proteinIdx = proteins.size();
            proteinIdMap.put(safeProteinId, proteinIdx);
            
            Isoform protein = new Isoform();
            protein.setId(safeProteinId);
            protein.setName(safeProteinId);
            protein.setGene(geneName);
            protein.setDescription("");
            protein.setSwissprot(true);
            proteins.add(protein);
            
            names.add(safeProteinId);
        } else {
            proteinIdx = proteinIdMap.get(safeProteinId);
        }
        
        if (!geneNameMap.containsKey(geneName)) {
            geneNameMap.put(geneName, genes.size());
            genes.add(geneName);
        }
        
        proteinToPrecursors.computeIfAbsent(proteinIdx, k -> new ArrayList<>()).add(precursorIndex);
        
        String precursorId = modifiedPeptideSequence + "/" + precursorCharge;
        precursors.add(precursorId);
        
        LibraryEntry libEntry = new LibraryEntry();
        libEntry.setName(precursorId);
        libEntry.setPidIndex(proteinIdx);
        libEntry.setProteotypic(proteotypic);
        
        Precursor precursor = new Precursor();
        precursor.setIndex(precursorIndex);
        precursor.setCharge(precursorCharge);
        precursor.setLength(peptideSequence != null ? peptideSequence.length() : 0);
        precursor.setMz(precursorMz);
        precursor.setiRT(normalizedRetentionTime);
        precursor.setiIM(precursorIonMobility);
        
        precursor.setFragments(fragments);
        
        libEntry.setTarget(precursor);
        entries.add(libEntry);
    }
    
    public void convertAndWrite(String outputPath) throws IOException {
        SpectralLibrary library = convert();
        
        long startWrite = System.nanoTime();
        DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
        writer.write(outputPath);
        long endWrite = System.nanoTime();
        System.out.println("Write time: " + (endWrite - startWrite) / 1_000_000 + " ms");
    }

    public static void main(String[] args) {
        String parquetFilePath = "I:\\test_transfer_learning\\dia_1\\fragpipe-predicted-speclib.parquet";
        String outputFilePath = "I:\\test_transfer_learning\\dia_1\\fragpipe-predicted-speclib.speclib";

        long startTime = System.nanoTime();
        ParquetToSpecLib converter = new ParquetToSpecLib(parquetFilePath);
        try {
            long convertStart = System.nanoTime();
            SpectralLibrary library = converter.convert();
            long convertEnd = System.nanoTime();
            System.out.println("Conversion time: " + (convertEnd - convertStart) / 1_000_000 + " ms");
            
            long writeStart = System.nanoTime();
            DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
            writer.write(outputFilePath);
            long writeEnd = System.nanoTime();
            System.out.println("Final write time: " + (writeEnd - writeStart) / 1_000_000 + " ms");
            
            System.out.println("Library stats: " + library.getEntries().size() + " precursors, " + 
                             library.getProteins().size() + " proteins");
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
