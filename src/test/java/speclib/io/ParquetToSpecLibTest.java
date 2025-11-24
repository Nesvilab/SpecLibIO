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

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.sql.*;
import java.util.*;

import static org.junit.Assert.*;

public class ParquetToSpecLibTest {

    @SuppressWarnings("unused")
    private static class ParquetRow {
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

    private File tempSpecLibFile;

    @Before
    public void setUp() throws IOException {
        tempSpecLibFile = File.createTempFile("test-speclib-", ".speclib");
    }

    @After
    public void tearDown() {
        if (tempSpecLibFile != null && tempSpecLibFile.exists()) {
            tempSpecLibFile.delete();
        }
    }

    @Test
    public void testRoundTripConversion() throws Exception {
        String parquetPath = getResourcePath("test.parquet");
        
        Map<String, List<ParquetRow>> originalData = readParquetData(parquetPath);
        
        Map<String, String> proteinToGeneMap = new HashMap<>();
        for (List<ParquetRow> rows : originalData.values()) {
            if (!rows.isEmpty() && rows.get(0).proteinId != null) {
                String proteinName = rows.get(0).proteinId;
                proteinToGeneMap.put(proteinName, "GENE_" + proteinName);
            }
        }
        
        ParquetToSpecLib converter = new ParquetToSpecLib(parquetPath, proteinToGeneMap, -1, true, true, true);
        SpectralLibrary library = converter.convert();
        
        verifyLibraryAgainstParquet(library, originalData, proteinToGeneMap);
        
        for (Isoform protein : library.getProteins()) {
            String expectedGene = "GENE_" + protein.getId();
            if (proteinToGeneMap.containsKey(protein.getId())) {
                assertEquals("Gene name should be preserved in round trip", expectedGene, protein.getGene());
            }
        }
    }

    @Test
    public void testPrecursorGrouping() throws Exception {
        String parquetPath = getResourcePath("test.parquet");
        
        Map<String, List<ParquetRow>> originalData = readParquetData(parquetPath);
        
        ParquetToSpecLib converter = new ParquetToSpecLib(parquetPath, Collections.emptyMap(), -1, true, true, true);
        SpectralLibrary library = converter.convert();
        
        Set<String> expectedPrecursors = originalData.keySet();
        assertEquals("Number of precursors should match", 
                     expectedPrecursors.size(), 
                     library.getEntries().size());
    }

    @Test
    public void testFragmentData() throws Exception {
        String parquetPath = getResourcePath("test.parquet");
        
        Map<String, List<ParquetRow>> originalData = readParquetData(parquetPath);
        
        ParquetToSpecLib converter = new ParquetToSpecLib(parquetPath, Collections.emptyMap(), -1, true, true, true);
        SpectralLibrary library = converter.convert();
        
        for (LibraryEntry entry : library.getEntries()) {
            Precursor precursor = entry.getTarget();
            String entryName = entry.getName();

            List<ParquetRow> expectedFragments = originalData.get(entryName);
            assertNotNull("Should have original data for precursor: " + entryName, expectedFragments);
            
            List<Product> actualFragments = precursor.getFragments();
            assertEquals("Fragment count should match for " + entryName,
                        expectedFragments.size(),
                        actualFragments.size());
        }
    }

    @Test
    public void testRetentionTimeRange() throws Exception {
        String parquetPath = getResourcePath("test.parquet");
        
        Map<String, List<ParquetRow>> originalData = readParquetData(parquetPath);
        
        float minRT = Float.MAX_VALUE;
        float maxRT = Float.MIN_VALUE;
        
        for (List<ParquetRow> rows : originalData.values()) {
            if (!rows.isEmpty()) {
                float rt = rows.get(0).normalizedRetentionTime;
                minRT = Math.min(minRT, rt);
                maxRT = Math.max(maxRT, rt);
            }
        }
        
        ParquetToSpecLib converter = new ParquetToSpecLib(parquetPath, Collections.emptyMap(), -1, true, true, true);
        SpectralLibrary library = converter.convert();
        
        assertEquals("Min RT should match", minRT, library.getiRTMin(), 0.0001);
        assertEquals("Max RT should match", maxRT, library.getiRTMax(), 0.0001);
    }

    @Test
    public void testProteinMapping() throws Exception {
        String parquetPath = getResourcePath("test.parquet");
        
        Map<String, List<ParquetRow>> originalData = readParquetData(parquetPath);
        
        Set<String> expectedProteins = new HashSet<>();
        Map<String, String> proteinToGeneMap = new HashMap<>();
        
        for (List<ParquetRow> rows : originalData.values()) {
            if (!rows.isEmpty() && rows.get(0).proteinId != null) {
                String proteinName = rows.get(0).proteinId;
                expectedProteins.add(proteinName);
                proteinToGeneMap.put(proteinName, "GENE_" + proteinName);
            }
        }
        
        ParquetToSpecLib converter = new ParquetToSpecLib(parquetPath, proteinToGeneMap, -1, true, true, true);
        SpectralLibrary library = converter.convert();
        
        assertEquals("Number of proteins should match",
                    expectedProteins.size(),
                    library.getProteins().size());
                    
        for (Isoform protein : library.getProteins()) {
            if (proteinToGeneMap.containsKey(protein.getId())) {
                assertEquals("Gene name should match map", 
                           proteinToGeneMap.get(protein.getId()), 
                           protein.getGene());
            } else {
                assertEquals("Gene name should default to protein ID", 
                           protein.getId(), 
                           protein.getGene());
            }
        }
    }

    private void verifyLibraryAgainstParquet(SpectralLibrary library, Map<String, List<ParquetRow>> originalData, Map<String, String> proteinToGeneMap) {
        assertNotNull("Library should not be null", library);
        assertNotNull("Library entries should not be null", library.getEntries());
        
        for (LibraryEntry entry : library.getEntries()) {
            Precursor precursor = entry.getTarget();
            String entryName = entry.getName();
            
            List<ParquetRow> originalRows = originalData.get(entryName);
            assertNotNull("Should have original data for precursor: " + entryName, originalRows);
            assertTrue("Should have at least one fragment", !originalRows.isEmpty());
            
            ParquetRow firstRow = originalRows.get(0);
            
            assertEquals("Precursor charge should match", 
                        firstRow.precursorCharge, 
                        precursor.getCharge());
            
            assertEquals("Precursor m/z should match", 
                        firstRow.precursorMz, 
                        precursor.getMz(), 
                        0.0001f);
            
            assertEquals("Retention time should match", 
                        firstRow.normalizedRetentionTime, 
                        precursor.getiRT(), 
                        0.0001f);
            
            float expectedIonMobility = parseIonMobility(firstRow.precursorIonMobility);
            assertEquals("Ion mobility should match", 
                        expectedIonMobility, 
                        precursor.getiIM(), 
                        0.0001f);
            
            assertEquals("Peptide length should match", 
                        firstRow.peptideSequence != null ? firstRow.peptideSequence.length() : 0, 
                        precursor.getLength());

            assertEquals("Modified peptide sequence should match",
                        firstRow.modifiedPeptideSequence,
                        entryName.substring(0, entryName.length() - 1));
            
            Isoform protein = library.getProteins().get(entry.getPidIndex());
            String expectedProteinId = firstRow.proteinId != null ? firstRow.proteinId : "";
            assertEquals("Protein ID should match",
                        expectedProteinId,
                        protein.getId());

            String expectedGene = proteinToGeneMap != null ? 
                proteinToGeneMap.getOrDefault(expectedProteinId, expectedProteinId) : 
                expectedProteinId;
            assertEquals("Gene name should match",
                        expectedGene,
                        protein.getGene());
            
            assertEquals("Proteotypic flag should match",
                        firstRow.proteotypic,
                        entry.getProteotypic());
            
            List<Product> fragments = precursor.getFragments();
            assertEquals("Fragment count should match", 
                        originalRows.size(), 
                        fragments.size());
            
            for (int i = 0; i < fragments.size(); i++) {
                Product fragment = fragments.get(i);
                ParquetRow expectedRow = originalRows.get(i);
                
                assertEquals("Fragment m/z should match at index " + i, 
                            expectedRow.productMz, 
                            fragment.getMz(), 
                            0.0001f);
                
                assertEquals("Fragment intensity should match at index " + i, 
                            expectedRow.libraryIntensity, 
                            fragment.getHeight(), 
                            0.0001f);
                
                assertEquals("Fragment charge should match at index " + i, 
                            (byte) expectedRow.fragmentCharge, 
                            fragment.getCharge());
                
                assertEquals("Fragment series number should match at index " + i,
                            (byte) expectedRow.fragmentSeriesNumber,
                            fragment.getIndex());
                
                byte expectedType = parseFragmentType(expectedRow.fragmentType);
                assertEquals("Fragment type should match at index " + i,
                            expectedType,
                            fragment.getType());
                
                byte expectedLoss = parseLossType(expectedRow.fragmentLossType);
                assertEquals("Fragment loss type should match at index " + i,
                            expectedLoss,
                            fragment.getLoss());
            }
        }
    }

    private Map<String, List<ParquetRow>> readParquetData(String parquetPath) throws SQLException {
        Map<String, List<ParquetRow>> data = new LinkedHashMap<>();
        
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:")) {
            String query = "SELECT " +
                "PrecursorMz, ProductMz, Annotation, ProteinId, " +
                "PeptideSequence, ModifiedPeptideSequence, PrecursorCharge, " +
                "LibraryIntensity, NormalizedRetentionTime, PrecursorIonMobility, " +
                "FragmentType, FragmentCharge, FragmentSeriesNumber, " +
                "FragmentLossType, Proteotypic " +
                "FROM read_parquet(?)";
            
            try (PreparedStatement stmt = conn.prepareStatement(query)) {
                stmt.setString(1, parquetPath);
                
                try (ResultSet rs = stmt.executeQuery()) {
                    while (rs.next()) {
                        ParquetRow row = new ParquetRow();
                        row.precursorMz = rs.getFloat("PrecursorMz");
                        row.productMz = rs.getFloat("ProductMz");
                        row.annotation = rs.getString("Annotation");
                        row.proteinId = rs.getString("ProteinId");
                        row.peptideSequence = rs.getString("PeptideSequence");
                        row.modifiedPeptideSequence = rs.getString("ModifiedPeptideSequence");
                        row.precursorCharge = rs.getShort("PrecursorCharge");
                        row.libraryIntensity = rs.getFloat("LibraryIntensity");
                        row.normalizedRetentionTime = rs.getFloat("NormalizedRetentionTime");
                        row.precursorIonMobility = rs.getString("PrecursorIonMobility");
                        row.fragmentType = rs.getString("FragmentType");
                        row.fragmentCharge = rs.getShort("FragmentCharge");
                        row.fragmentSeriesNumber = rs.getShort("FragmentSeriesNumber");
                        row.fragmentLossType = rs.getString("FragmentLossType");
                        row.proteotypic = rs.getShort("Proteotypic");
                        
                        String precursorKey = row.modifiedPeptideSequence + row.precursorCharge;
                        data.computeIfAbsent(precursorKey, k -> new ArrayList<>()).add(row);
                    }
                }
            }
        }
        
        return data;
    }

    private String getResourcePath(String resourceName) {
        ClassLoader classLoader = getClass().getClassLoader();
        File file = new File(classLoader.getResource(resourceName).getFile());
        return file.getAbsolutePath();
    }

    private static float parseIonMobility(String ionMobility) {
        if (ionMobility == null || ionMobility.isEmpty()) {
            return 0.0f;
        }
        try {
            return Float.parseFloat(ionMobility);
        } catch (NumberFormatException e) {
            return 0.0f;
        }
    }

    private static byte parseFragmentType(String fragmentType) {
        if (fragmentType == null || fragmentType.isEmpty()) {
            throw new IllegalArgumentException("Fragment type cannot be null or empty");
        }
        
        char firstChar = fragmentType.charAt(0);
        if (firstChar == 'b' || firstChar == 'B') {
            return 1;
        } else if (firstChar == 'y' || firstChar == 'Y') {
            return 2;
        }
        
        throw new IllegalArgumentException(
            "Unsupported fragment type: '" + fragmentType + "'. " +
            "Only 'b' and 'y' ion types are supported by DIA-NN speclib format.");
    }

    private static byte parseLossType(String lossType) {
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
}

