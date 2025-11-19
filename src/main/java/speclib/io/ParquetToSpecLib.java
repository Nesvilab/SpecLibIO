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
    
    @SuppressWarnings("unused")
    private static class FragmentRow {
        float precursorMz;
        float productMz;
        String annotation;
        String proteinId;
        String peptideSequence;
        String modifiedPeptideSequence;
        short precursorCharge;
        float libraryIntensity;
        float normalizedRetentionTime;
        String precursorIonMobility;
        String fragmentType;
        short fragmentCharge;
        short fragmentSeriesNumber;
        String fragmentLossType;
        short proteotypic;
    }
    
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
        List<FragmentRow> rows = readParquetFile();
        return buildSpectralLibrary(rows);
    }
    
    private List<FragmentRow> readParquetFile() throws IOException {
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:")) {            
            long rowCount = getRowCount(conn);
            List<FragmentRow> rows = new ArrayList<>((int) Math.min(rowCount, Integer.MAX_VALUE));
            
            String query = "SELECT " +
                "PrecursorMz, ProductMz, Annotation, ProteinId, " +
                "PeptideSequence, ModifiedPeptideSequence, PrecursorCharge, " +
                "LibraryIntensity, NormalizedRetentionTime, PrecursorIonMobility, " +
                "FragmentType, FragmentCharge, FragmentSeriesNumber, " +
                "FragmentLossType, Proteotypic " +
                "FROM read_parquet(?)";
            
            try (PreparedStatement stmt = conn.prepareStatement(query,
                    ResultSet.TYPE_FORWARD_ONLY, 
                    ResultSet.CONCUR_READ_ONLY)) {
                stmt.setFetchSize(1 << 10);
                stmt.setString(1, parquetFilePath);
                
                try (ResultSet rs = stmt.executeQuery()) {
                    while (rs.next()) {
                        FragmentRow row = new FragmentRow();
                        row.precursorMz = rs.getFloat(1);
                        row.productMz = rs.getFloat(2);
                        row.annotation = rs.getString(3);
                        row.proteinId = rs.getString(4);
                        row.peptideSequence = rs.getString(5);
                        row.modifiedPeptideSequence = rs.getString(6);
                        row.precursorCharge = rs.getShort(7);
                        row.libraryIntensity = rs.getFloat(8);
                        row.normalizedRetentionTime = rs.getFloat(9);
                        row.precursorIonMobility = rs.getString(10);
                        row.fragmentType = rs.getString(11);
                        row.fragmentCharge = rs.getShort(12);
                        row.fragmentSeriesNumber = rs.getShort(13);
                        row.fragmentLossType = rs.getString(14);
                        row.proteotypic = rs.getShort(15);
                        rows.add(row);
                    }
                }
            }
            
            return rows;
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

    private SpectralLibrary buildSpectralLibrary(List<FragmentRow> rows) {
        SpectralLibrary library = new SpectralLibrary();
        
        int estimatedPrecursors = rows.size() / 16;
        int estimatedProteins = Math.max(1 << 13, rows.size() / 128);
        
        Map<String, Integer> proteinIdMap = new LinkedHashMap<>(estimatedProteins);
        Map<String, Integer> geneNameMap = new LinkedHashMap<>(estimatedProteins);
        Map<String, Integer> precursorMap = new LinkedHashMap<>(estimatedPrecursors);
        
        List<Isoform> proteins = new ArrayList<>(estimatedProteins);
        List<ProteinGroup> proteinIds = new ArrayList<>(estimatedProteins);
        List<String> precursors = new ArrayList<>(estimatedPrecursors);
        List<String> names = new ArrayList<>(estimatedProteins);
        List<String> genes = new ArrayList<>(estimatedProteins);
        
        Map<String, List<FragmentRow>> precursorFragments = new LinkedHashMap<>(estimatedPrecursors);
        for (FragmentRow row : rows) {
            String precursorKey = row.modifiedPeptideSequence + "_" + row.precursorCharge;
            precursorFragments.computeIfAbsent(precursorKey, k -> new ArrayList<>(16)).add(row);
        }
        
        List<LibraryEntry> entries = new ArrayList<>(estimatedPrecursors);
        
        double minRT = Double.MAX_VALUE;
        double maxRT = Double.MIN_VALUE;
        
        int precursorIndex = 0;
        for (Map.Entry<String, List<FragmentRow>> entry : precursorFragments.entrySet()) {
            List<FragmentRow> fragments = entry.getValue();
            if (fragments.isEmpty()) {
                continue;
            }
            
            FragmentRow firstRow = fragments.get(0);
            
            String proteinId = firstRow.proteinId != null ? firstRow.proteinId : "";
            String geneName = proteinToGeneMap.getOrDefault(proteinId, proteinId);
            
            int proteinIdx;
            if (!proteinIdMap.containsKey(proteinId)) {
                proteinIdx = proteins.size();
                proteinIdMap.put(proteinId, proteinIdx);
                
                Isoform protein = new Isoform();
                protein.setId(proteinId);
                protein.setName(proteinId);
                protein.setGene(geneName);
                protein.setDescription("");
                protein.setSwissprot(true);
                proteins.add(protein);
            } else {
                proteinIdx = proteinIdMap.get(proteinId);
            }
            
            if (!geneNameMap.containsKey(geneName)) {
                geneNameMap.put(geneName, genes.size());
                genes.add(geneName);
            }
            
            if (!proteinIdMap.containsKey(proteinId)) {
                names.add(proteinId);
            }
            
            precursorMap.put(entry.getKey(), precursorIndex);
            precursors.add(firstRow.modifiedPeptideSequence);
            
            LibraryEntry libEntry = new LibraryEntry();
            libEntry.setName(firstRow.modifiedPeptideSequence);
            libEntry.setPidIndex(proteinIdx);
            libEntry.setProteotypic(firstRow.proteotypic);
            
            Peptide peptide = new Peptide();
            peptide.setIndex(precursorIndex);
            peptide.setCharge(firstRow.precursorCharge);
            peptide.setLength(firstRow.peptideSequence != null ? firstRow.peptideSequence.length() : 0);
            peptide.setMz(firstRow.precursorMz);
            peptide.setiRT(firstRow.normalizedRetentionTime);
            float ionMobility = parseIonMobility(firstRow.precursorIonMobility);
            peptide.setiIM(ionMobility);
            
            if (firstRow.normalizedRetentionTime < minRT) {
                minRT = firstRow.normalizedRetentionTime;
            }
            if (firstRow.normalizedRetentionTime > maxRT) {
                maxRT = firstRow.normalizedRetentionTime;
            }
            
            List<Product> products = new ArrayList<>();
            for (FragmentRow fragRow : fragments) {
                Product product = new Product();
                product.setMz(fragRow.productMz);
                product.setHeight(fragRow.libraryIntensity);
                product.setCharge((byte) fragRow.fragmentCharge);
                product.setType(parseFragmentType(fragRow.fragmentType));
                product.setIndex((byte) fragRow.fragmentSeriesNumber);
                product.setLoss(parseLossType(fragRow.fragmentLossType));
                products.add(product);
            }
            peptide.setFragments(products);
            
            libEntry.setTarget(peptide);
            entries.add(libEntry);
            
            precursorIndex++;
        }
        
        for (int i = 0; i < proteins.size(); i++) {
            Isoform protein = proteins.get(i);
            String geneName = protein.getGene();
            
            protein.setNameIndex(i < names.size() ? i : 0);
            protein.setGeneIndex(geneNameMap.getOrDefault(geneName, 0));
            
            Set<Integer> precursorIndices = new LinkedHashSet<>();
            for (int j = 0; j < entries.size(); j++) {
                if (entries.get(j).getPidIndex() == i) {
                    precursorIndices.add(j);
                }
            }
            protein.setPrecursors(precursorIndices);
        }
        
        for (Isoform protein : proteins) {
            ProteinGroup pg = new ProteinGroup();
            pg.setIds(protein.getId());
            pg.setNames(protein.getName());
            pg.setGenes(protein.getGene());
            
            List<Integer> precursorIndices = new ArrayList<>(protein.getPrecursors());
            pg.setPrecursors(precursorIndices);
            
            Set<Integer> proteinSet = new LinkedHashSet<>();
            proteinSet.add(proteinIdMap.get(protein.getId()));
            pg.setProteins(proteinSet);
            
            List<Integer> nameIndices = new ArrayList<>();
            nameIndices.add(protein.getNameIndex());
            pg.setNameIndices(nameIndices);
            
            List<Integer> geneIndices = new ArrayList<>();
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
    
    public static byte parseFragmentType(String fragmentType) {
        if (fragmentType == null || fragmentType.isEmpty()) {
            throw new IllegalArgumentException("Fragment type cannot be null or empty");
        }
        
        switch (fragmentType.toLowerCase()) {
            case "b":
                return 1;
            case "y":
                return 2;
            default:
                throw new IllegalArgumentException(
                    "Unsupported fragment type: '" + fragmentType + "'. " +
                    "Only 'b' and 'y' ion types are supported by DIA-NN speclib format.");
        }
    }
    
    public static byte parseLossType(String lossType) {
        if (lossType == null || lossType.isEmpty()) {
            return 0;
        }
        
        switch (lossType.toLowerCase()) {
            case "noloss":
                return 0;
            case "h2o":
            case "-h2o":
                return 1;
            case "nh3":
            case "-nh3":
                return 2;
            case "co":
            case "-co":
                return 3;
            case "n":
            case "-n":
                return 4;
            case "other":
            case "-other":
                return 5;
            default:
                return 0;
        }
    }
    
    public static float parseIonMobility(String ionMobility) {
        if (ionMobility == null || ionMobility.isEmpty()) {
            return 0.0f;
        }
        try {
            return Float.parseFloat(ionMobility);
        } catch (NumberFormatException e) {
            return 0.0f;
        }
    }
    
    public void convertAndWrite(String outputPath) throws IOException {
        SpectralLibrary library = convert();
        DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
        writer.write(outputPath);
    }
}
