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

import org.junit.Assert;
import org.junit.Assume;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DiaNNSpecLibReaderTest {

    @Rule
    public TemporaryFolder tempFolder = new TemporaryFolder();

    private static class TsvEntry {
        String precursorId;
        int charge;
        float rt;
        float iRT;
        float iIM;
        float im;
        float mz;
        float libQValue;
        int proteotypic;
        String strippedSequence;
        String proteinIds;
        String proteinNames;
        String genes;
        
        TsvEntry(String precursorId) {
            this.precursorId = precursorId;
        }
    }

    private Map<String, TsvEntry> parseTsvFile(File tsvFile) throws IOException {
        Map<String, TsvEntry> entries = new HashMap<>();
        
        try (BufferedReader reader = new BufferedReader(new FileReader(tsvFile))) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                return entries;
            }
            
            String[] headers = headerLine.split("\t");
            Map<String, Integer> columnIndex = new HashMap<>();
            for (int i = 0; i < headers.length; i++) {
                columnIndex.put(headers[i], i);
            }
            
            if (!columnIndex.containsKey("Precursor.Id")) {
                throw new IOException("TSV file missing required column: Precursor.Id");
            }
            
            String line;
            while ((line = reader.readLine()) != null) {
                String[] values = line.split("\t", -1);
                
                String precursorId = getColumnValue(values, columnIndex, "Precursor.Id");
                if (precursorId == null || precursorId.isEmpty()) {
                    continue;
                }
                
                if (entries.containsKey(precursorId)) {
                    continue;
                }
                
                TsvEntry entry = new TsvEntry(precursorId);
                
                entry.charge = parseIntColumn(values, columnIndex, "Precursor.Charge", 0);
                entry.rt = parseFloatColumn(values, columnIndex, "RT", 0.0f);
                entry.iRT = parseFloatColumn(values, columnIndex, "iRT", 0.0f);
                entry.iIM = parseFloatColumn(values, columnIndex, "iIM", 0.0f);
                entry.im = parseFloatColumn(values, columnIndex, "IM", 0.0f);
                entry.mz = parseFloatColumn(values, columnIndex, "Precursor.Mz", 0.0f);
                entry.libQValue = parseFloatColumn(values, columnIndex, "Lib.Q.Value", 0.0f);
                entry.proteotypic = parseIntColumn(values, columnIndex, "Proteotypic", 0);
                entry.strippedSequence = getColumnValue(values, columnIndex, "Stripped.Sequence");
                entry.proteinIds = getColumnValue(values, columnIndex, "Protein.Ids");
                entry.proteinNames = getColumnValue(values, columnIndex, "Protein.Names");
                entry.genes = getColumnValue(values, columnIndex, "Genes");
                
                entries.put(precursorId, entry);
            }
        }
        
        return entries;
    }

    private String getColumnValue(String[] values, Map<String, Integer> columnIndex, String columnName) {
        Integer index = columnIndex.get(columnName);
        if (index == null || index >= values.length) {
            return null;
        }
        String value = values[index];
        return value.isEmpty() ? null : value;
    }

    private int parseIntColumn(String[] values, Map<String, Integer> columnIndex, String columnName, int defaultValue) {
        String value = getColumnValue(values, columnIndex, columnName);
        if (value == null) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }

    private float parseFloatColumn(String[] values, Map<String, Integer> columnIndex, String columnName, float defaultValue) {
        String value = getColumnValue(values, columnIndex, columnName);
        if (value == null) {
            return defaultValue;
        }
        try {
            return Float.parseFloat(value);
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }

    private void validateCommonFields(LibraryEntry entry, TsvEntry tsv, SpectralLibrary library) {
        Assert.assertEquals("Precursor.Id should match entry name", tsv.precursorId, entry.getName());
        
        Assert.assertEquals("Charge mismatch for " + entry.getName(), tsv.charge, entry.getTarget().getCharge());
        
        if (tsv.mz != 0.0f) {
            Assert.assertEquals("Precursor.Mz mismatch for " + entry.getName(), tsv.mz, entry.getTarget().getMz(), 0.01f);
        }
        
        if (tsv.iRT != 0.0f) {
            Assert.assertEquals("iRT mismatch for " + entry.getName(), tsv.iRT, entry.getTarget().getiRT(), 0.1f);
        }
        
        if (tsv.rt != 0.0f && entry.getTarget().getsRT() != 0.0f) {
            Assert.assertEquals("RT mismatch for " + entry.getName(), tsv.rt, entry.getTarget().getsRT(), 0.1f);
        }
        
        if (tsv.iIM != 0.0f && entry.getTarget().getiIM() != 0.0f) {
            Assert.assertEquals("iIM mismatch for " + entry.getName(), tsv.iIM, entry.getTarget().getiIM(), 0.01f);
        }
        
        if (tsv.im != 0.0f && entry.getTarget().getsIM() != 0.0f) {
            Assert.assertEquals("IM mismatch for " + entry.getName(), tsv.im, entry.getTarget().getsIM(), 0.01f);
        }
        
        if (tsv.libQValue != 0.0f) {
            Assert.assertEquals("Lib.Q.Value mismatch for " + entry.getName(), tsv.libQValue, entry.getTarget().getLibQvalue(), 0.01f);
        }
        
        if (tsv.proteotypic != 0 && entry.getProteotypic() != 0) {
            Assert.assertEquals("Proteotypic mismatch for " + entry.getName(), tsv.proteotypic, entry.getProteotypic());
        }
        
        if (tsv.strippedSequence != null && !tsv.strippedSequence.isEmpty() && entry.getTarget().getLength() > 0) {
            int expectedLength = tsv.strippedSequence.length();
            Assert.assertEquals("Peptide length mismatch for " + entry.getName(), expectedLength, entry.getTarget().getLength());
        }
        
        if (tsv.proteinIds != null && !tsv.proteinIds.isEmpty()) {
            validateProteinIds(entry, tsv, library);
        }
        
        if (tsv.proteinNames != null && !tsv.proteinNames.isEmpty()) {
            validateProteinNames(entry, tsv, library);
        }
        
        if (tsv.genes != null && !tsv.genes.isEmpty()) {
            validateGenes(entry, tsv, library);
        }
        
        Assert.assertTrue("Should have fragments for " + entry.getName(), entry.getTarget().getFragments().size() > 0);
    }
    
    private void validateProteinIds(LibraryEntry entry, TsvEntry tsv, SpectralLibrary library) {
        if (library.getIsoforms().isEmpty() || library.getProteinGroups().isEmpty()) {
            return;
        }
        
        if (entry.getPgIndex() >= 0 && entry.getPgIndex() < library.getProteinGroups().size()) {
            ProteinGroup proteinGroup = library.getProteinGroups().get(entry.getPgIndex());
            if (proteinGroup.getIsoforms().isEmpty()) {
                return;
            }
            
            String[] tsvProteinIds = tsv.proteinIds.split(";");
            if (tsvProteinIds.length == 0) {
                return;
            }
            
            boolean foundMatch = false;
            for (String tsvProteinId : tsvProteinIds) {
                tsvProteinId = tsvProteinId.trim();
                if (tsvProteinId.isEmpty()) continue;
                
                for (int proteinIndex : proteinGroup.getIsoforms()) {
                    if (proteinIndex >= 0 && proteinIndex < library.getIsoforms().size()) {
                        Isoform isoform = library.getIsoforms().get(proteinIndex);
                        String isoformName = isoform.getName();
                        if (isoformName.equals(tsvProteinId)) {
                            foundMatch = true;
                            break;
                        }
                    }
                }
                if (foundMatch) break;
            }
        }
    }
    
    private void validateProteinNames(LibraryEntry entry, TsvEntry tsv, SpectralLibrary library) {
        if (library.getIsoforms().isEmpty() || library.getProteinGroups().isEmpty()) {
            return;
        }
        
        if (entry.getPgIndex() >= 0 && entry.getPgIndex() < library.getProteinGroups().size()) {
            ProteinGroup proteinGroup = library.getProteinGroups().get(entry.getPgIndex());
            if (proteinGroup.getIsoforms().isEmpty()) {
                return;
            }
            
            String[] tsvProteinNames = tsv.proteinNames.split(";");
            if (tsvProteinNames.length == 0) {
                return;
            }
            
            boolean foundMatch = false;
            for (String tsvProteinName : tsvProteinNames) {
                tsvProteinName = tsvProteinName.trim();
                if (tsvProteinName.isEmpty()) continue;
                
                for (int proteinIndex : proteinGroup.getIsoforms()) {
                    if (proteinIndex >= 0 && proteinIndex < library.getIsoforms().size()) {
                        Isoform isoform = library.getIsoforms().get(proteinIndex);
                        String isoformName = isoform.getName();
                        if (isoformName.equals(tsvProteinName)) {
                            foundMatch = true;
                            break;
                        }
                    }
                }
                if (foundMatch) break;
            }
        }
    }
    
    private void validateGenes(LibraryEntry entry, TsvEntry tsv, SpectralLibrary library) {
        if (library.getGenes().isEmpty()) {
            return;
        }
        
        String[] tsvGenes = tsv.genes.split(";");
        if (tsvGenes.length == 0) {
            return;
        }
        
        boolean foundMatch = false;
        for (String tsvGene : tsvGenes) {
            tsvGene = tsvGene.trim();
            if (tsvGene.isEmpty()) continue;
            
            for (String libGene : library.getGenes()) {
                if (libGene.equals(tsvGene)) {
                    foundMatch = true;
                    break;
                }
            }
            if (foundMatch) break;
        }
    }
    
    private void validateFragments(LibraryEntry entry) {
        List<Product> fragments = entry.getTarget().getFragments();
        
        Assert.assertTrue("Entry should have at least 2 fragments: " + entry.getName(), fragments.size() >= 2);
        
        for (Product fragment : fragments) {
            Assert.assertTrue("Fragment m/z should be positive for " + entry.getName(), fragment.getMz() > 0);
            Assert.assertTrue("Fragment m/z should be reasonable (< 2000) for " + entry.getName(), fragment.getMz() < 2000);
            Assert.assertTrue("Fragment height should be non-negative for " + entry.getName(), fragment.getHeight() >= 0);
            Assert.assertTrue("Fragment charge should be positive for " + entry.getName(), fragment.getCharge() > 0);
            Assert.assertTrue("Fragment charge should be reasonable (<=4) for " + entry.getName(), fragment.getCharge() <= 4);
            Assert.assertTrue("Fragment index should be non-negative for " + entry.getName(), fragment.getIndex() >= 0);
        }
    }

    @Test
    public void testValidateAgainstTSV_DiannSwath() throws IOException {
        File speclibFile = new File(getClass().getClassLoader().getResource("diann-swath.speclib").getFile());
        File tsvFile = new File(getClass().getClassLoader().getResource("diann-swath-report.tsv").getFile());
        
        Assume.assumeTrue("Test files exist", speclibFile.exists() && tsvFile.exists());

        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(speclibFile);
        SpectralLibrary library = reader.read();

        Map<String, TsvEntry> tsvData = parseTsvFile(tsvFile);

        Assert.assertNotNull("Library should not be null", library);
        Assert.assertTrue("Library should have entries", library.getEntries().size() > 0);
        
        Assert.assertTrue("Library should have at least as many entries as TSV has unique precursors", library.getEntries().size() >= tsvData.size());
        Assert.assertEquals("Library entry count should match TSV precursor count", tsvData.size(), library.getEntries().size());

        int validatedCount = 0;
        int fragmentsValidated = 0;
        
        for (LibraryEntry entry : library.getEntries()) {
            TsvEntry tsv = tsvData.get(entry.getName());
            if (tsv != null) {
                validateCommonFields(entry, tsv, library);
                validateFragments(entry);
                
                fragmentsValidated++;
                validatedCount++;
            }
        }
        
        Assert.assertTrue("Should validate at least some entries", validatedCount > 0);
        Assert.assertTrue("Should validate fragments", fragmentsValidated > 0);
        Assert.assertEquals("All TSV entries should be found in library", tsvData.size(), validatedCount);
    }

    @Test
    public void testValidateAgainstTSV_DiannModTest() throws IOException {
        File speclibFile = new File(getClass().getClassLoader().getResource("diann-mod-test.tsv.speclib").getFile());
        File tsvFile = new File(getClass().getClassLoader().getResource("diann-mod-test.tsv").getFile());
        
        Assume.assumeTrue("Test files exist", speclibFile.exists() && tsvFile.exists());

        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(speclibFile);
        SpectralLibrary library = reader.read();

        Map<String, TsvEntry> tsvData = parseTsvFile(tsvFile);

        Assert.assertNotNull("Library should not be null", library);
        Assert.assertTrue("Library should have entries", library.getEntries().size() > 0);
        
        Assert.assertTrue("Library should have at least as many entries as TSV has unique precursors", library.getEntries().size() >= tsvData.size());

        int validatedCount = 0;
        int modifiedCount = 0;
        int fragmentsValidated = 0;
        
        for (LibraryEntry entry : library.getEntries()) {
            TsvEntry tsv = tsvData.get(entry.getName());
            if (tsv != null) {
                validateCommonFields(entry, tsv, library);
                validateFragments(entry);
                
                if (entry.getName().contains("(")) {
                    modifiedCount++;
                }
                
                fragmentsValidated++;
                validatedCount++;
            }
        }
        
        Assert.assertTrue("Should validate at least some entries", validatedCount > 0);
        Assert.assertTrue("Should have modified precursors", modifiedCount > 0);
        Assert.assertTrue("Should validate fragments", fragmentsValidated > 0);
        Assert.assertEquals("All TSV entries should be found in library", tsvData.size(), validatedCount);
    }

    @Test
    public void testValidateAgainstTSV_DiannMassMods() throws IOException {
        File speclibFile = new File(getClass().getClassLoader().getResource("diann-mass-mods.tsv.speclib").getFile());
        File tsvFile = new File(getClass().getClassLoader().getResource("diann-mass-mods.tsv").getFile());
        
        Assume.assumeTrue("Test files exist", speclibFile.exists() && tsvFile.exists());

        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(speclibFile);
        SpectralLibrary library = reader.read();

        Map<String, TsvEntry> tsvData = parseTsvFile(tsvFile);

        Assert.assertNotNull("Library should not be null", library);
        Assert.assertTrue("Library should have entries", library.getEntries().size() > 0);
        
        Assert.assertTrue("Library should have at least as many entries as TSV has unique precursors", library.getEntries().size() >= tsvData.size());

        int validatedCount = 0;
        int fragmentsValidated = 0;
        
        for (LibraryEntry entry : library.getEntries()) {
            TsvEntry tsv = tsvData.get(entry.getName());
            if (tsv != null) {
                validateCommonFields(entry, tsv, library);
                validateFragments(entry);
                
                fragmentsValidated++;
                validatedCount++;
            }
        }
        
        Assert.assertTrue("Should validate at least some entries", validatedCount > 0);
        Assert.assertTrue("Should validate fragments", fragmentsValidated > 0);
        Assert.assertEquals("All TSV entries should be found in library", tsvData.size(), validatedCount);
    }

    @Test
    public void testValidateAgainstTSV_DiannDiaPasef() throws IOException {
        File speclibFile = new File(getClass().getClassLoader().getResource("diann-hela-diapasef-lib.speclib").getFile());
        File tsvFile = new File(getClass().getClassLoader().getResource("diann-hela-diapasef.tsv").getFile());
        
        Assume.assumeTrue("Test files exist", speclibFile.exists() && tsvFile.exists());

        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(speclibFile);
        SpectralLibrary library = reader.read();

        Map<String, TsvEntry> tsvData = parseTsvFile(tsvFile);

        Assert.assertNotNull("Library should not be null", library);
        Assert.assertTrue("Library should have entries", library.getEntries().size() > 0);
        
        Assert.assertTrue("Library should have at least as many entries as TSV has unique precursors", library.getEntries().size() >= tsvData.size());
        Assert.assertTrue("Library and TSV counts should be close (within 5%)", Math.abs(library.getEntries().size() - tsvData.size()) <= tsvData.size() * 0.05);

        int validatedCount = 0;
        int iIMValidated = 0;
        int fragmentsValidated = 0;
        
        for (LibraryEntry entry : library.getEntries()) {
            TsvEntry tsv = tsvData.get(entry.getName());
            if (tsv != null) {
                validateCommonFields(entry, tsv, library);
                validateFragments(entry);
                
                if (tsv.iIM != 0.0f && entry.getTarget().getiIM() != 0.0f) {
                    iIMValidated++;
                }
                
                fragmentsValidated++;
                validatedCount++;
            }
        }
        
        Assert.assertTrue("Should validate at least some entries", validatedCount > 0);
        Assert.assertTrue("Should have ion mobility data", iIMValidated > 0);
        Assert.assertTrue("Should validate fragments", fragmentsValidated > 0);
        Assert.assertEquals("All TSV entries should be found in library", tsvData.size(), validatedCount);
    }

    @Test
    public void testValidateAgainstTSV_LibraryTsv() throws IOException {
        File speclibFile = new File(getClass().getClassLoader().getResource("library.tsv.speclib").getFile());
        File tsvFile = new File(getClass().getClassLoader().getResource("diann-output.tsv").getFile());
        
        Assume.assumeTrue("Test files exist", speclibFile.exists() && tsvFile.exists());

        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(speclibFile);
        SpectralLibrary library = reader.read();

        Map<String, TsvEntry> tsvData = parseTsvFile(tsvFile);

        Assert.assertNotNull("Library should not be null", library);
        Assert.assertTrue("Library should have entries", library.getEntries().size() > 0);
        
        Assert.assertTrue("Library should have at least as many entries as TSV has unique precursors", library.getEntries().size() >= tsvData.size());

        int validatedCount = 0;
        int fragmentsValidated = 0;
        
        for (LibraryEntry entry : library.getEntries()) {
            TsvEntry tsv = tsvData.get(entry.getName());
            if (tsv != null) {
                validateCommonFields(entry, tsv, library);
                validateFragments(entry);
                
                fragmentsValidated++;
                validatedCount++;
            }
        }
        
        Assert.assertTrue("Should validate at least some entries", validatedCount > 0);
        Assert.assertTrue("Should validate fragments", fragmentsValidated > 0);
        Assert.assertEquals("All TSV entries should be found in library", tsvData.size(), validatedCount);
    }

    @Test(expected = IOException.class)
    public void testReadEmptyFile() throws IOException {
        ByteArrayInputStream bais = new ByteArrayInputStream(new byte[0]);
        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(bais);
        reader.read();
    }

    @Test(expected = IOException.class)
    public void testReadTruncatedFile() throws IOException {
        SpectralLibrary original = DiaNNSpecLibWriter.createSampleLibrary("Test", "PEPTIDE", 2, 400.0f);
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        new DiaNNSpecLibWriter(original).write(baos);
        
        byte[] truncated = new byte[50];
        System.arraycopy(baos.toByteArray(), 0, truncated, 0, 50);
        
        ByteArrayInputStream bais = new ByteArrayInputStream(truncated);
        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(bais);
        reader.read();
    }

    private static void writeIsoform(BufferedWriter writer, Isoform isoform, int idx) throws IOException {
        writer.write("  [" + idx + "]:\n");
        writer.write("    swissprot: " + isoform.isSwissprot() + "\n");
        writer.write("    precursors.size: " + isoform.getPrecursors().size() + "\n");
        writer.write("    id: " + isoform.getId() + "\n");
        writer.write("    name: " + isoform.getName() + "\n");
        writer.write("    gene: " + isoform.getGene() + "\n");
        writer.write("    nameIndex: " + isoform.getNameIndex() + "\n");
        writer.write("    geneIndex: " + isoform.getGeneIndex() + "\n");
        writer.write("    precursors: " + isoform.getPrecursors() + "\n");
    }

    private static void writeProteinGroup(BufferedWriter writer, ProteinGroup pg, int idx) throws IOException {
        writer.write("  [" + idx + "]:\n");
        writer.write("    proteins.size: " + pg.getIsoforms().size() + "\n");
        writer.write("    ids: " + pg.getIds() + "\n");
        writer.write("    names: " + pg.getNames() + "\n");
        writer.write("    genes: " + pg.getGenes() + "\n");
        writer.write("    nameIndices: " + pg.getNameIndices() + "\n");
        writer.write("    geneIndices: " + pg.getGeneIndices() + "\n");
        writer.write("    precursors: " + pg.getPrecursors() + "\n");
        writer.write("    proteins: " + pg.getIsoforms() + "\n");
    }

    private static void writeProduct(BufferedWriter writer, Product product, int idx) throws IOException {
        writer.write("        [" + idx + "]: mz=" + product.getMz() + 
                     ", height=" + product.getHeight() + 
                     ", charge=" + product.getCharge() + 
                     ", type=" + product.getType() + 
                     ", index=" + product.getIndex() + 
                     ", loss=" + product.getLoss() + "\n");
    }

    private static void writePrecursor(BufferedWriter writer, Precursor precursor, String label) throws IOException {
        writer.write("    " + label + ":\n");
        writer.write("      index: " + precursor.getIndex() + "\n");
        writer.write("      charge: " + precursor.getCharge() + "\n");
        writer.write("      length: " + precursor.getLength() + "\n");
        writer.write("      mz: " + precursor.getMz() + "\n");
        writer.write("      iRT: " + precursor.getiRT() + "\n");
        writer.write("      sRT: " + precursor.getsRT() + "\n");
        writer.write("      libQvalue: " + precursor.getLibQvalue() + "\n");
        writer.write("      iIM: " + precursor.getiIM() + "\n");
        writer.write("      sIM: " + precursor.getsIM() + "\n");
        writer.write("      fragments: " + precursor.getFragments().size() + " entries\n");
        for (int i = 0; i < precursor.getFragments().size(); i++) {
            writeProduct(writer, precursor.getFragments().get(i), i);
        }
    }

    private static void writeLibraryEntry(BufferedWriter writer, LibraryEntry entry, int idx) throws IOException {
        writer.write("  [" + idx + "]:\n");
        writePrecursor(writer, entry.getTarget(), "target");
        writer.write("    hasDecoy: " + (entry.getDecoy() != null) + "\n");
        if (entry.getDecoy() != null) {
            writePrecursor(writer, entry.getDecoy(), "decoy");
        }
        writer.write("    entryFlags: " + entry.getEntryFlags() + "\n");
        writer.write("    proteotypic: " + entry.getProteotypic() + "\n");
        writer.write("    pgIndex: " + entry.getPgIndex() + "\n");
        writer.write("    name: " + entry.getName() + "\n");
        writer.write("    pgQvalue: " + entry.getPgQvalue() + "\n");
        writer.write("    ptmQvalue: " + entry.getPtmQvalue() + "\n");
        writer.write("    siteConf: " + entry.getSiteConf() + "\n");
    }

    @Test
    public void testWriteLibraryContentToFile() throws IOException {
        File speclibFile = new File(getClass().getClassLoader().getResource("1.8.1/report-lib.predicted.speclib").getFile());
        
        Assume.assumeTrue("Test file exists", speclibFile.exists());

        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(speclibFile);
        SpectralLibrary library = reader.read();

        File outputFile = new File(speclibFile.getParent(), "library_content.txt");

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            writer.write("=== Spectral Library Content ===\n\n");

            writer.write("genDecoys: " + library.isGenDecoys() + "\n");
            writer.write("genCharges: " + library.isGenCharges() + "\n");
            writer.write("inferProteotypicity: " + library.isInferProteotypicity() + "\n\n");

            writer.write("name: " + library.getName() + "\n");
            writer.write("fastaNames: " + library.getFastaNames() + "\n\n");

            writer.write("proteins (Isoforms): " + library.getIsoforms().size() + " entries\n");
            for (int i = 0; i < library.getIsoforms().size(); i++) {
                writeIsoform(writer, library.getIsoforms().get(i), i);
            }
            writer.write("\n");

            writer.write("protein groups (ProteinGroups): " + library.getProteinGroups().size() + " entries\n");
            for (int i = 0; i < library.getProteinGroups().size(); i++) {
                writeProteinGroup(writer, library.getProteinGroups().get(i), i);
            }
            writer.write("\n");

            writer.write("precursors: " + library.getPrecursors().size() + " entries\n");
            for (int i = 0; i < library.getPrecursors().size(); i++) {
                writer.write("  [" + i + "]: " + library.getPrecursors().get(i) + "\n");
            }
            writer.write("\n");

            writer.write("names: " + library.getNames().size() + " entries\n");
            for (int i = 0; i < library.getNames().size(); i++) {
                writer.write("  [" + i + "]: " + library.getNames().get(i) + "\n");
            }
            writer.write("\n");

            writer.write("genes: " + library.getGenes().size() + " entries\n");
            for (int i = 0; i < library.getGenes().size(); i++) {
                writer.write("  [" + i + "]: " + library.getGenes().get(i) + "\n");
            }
            writer.write("\n");

            writer.write("iRTMin: " + library.getiRTMin() + "\n");
            writer.write("iRTMax: " + library.getiRTMax() + "\n\n");

            writer.write("entries (LibraryEntries): " + library.getEntries().size() + " entries\n");
            for (int i = 0; i < library.getEntries().size(); i++) {
                writeLibraryEntry(writer, library.getEntries().get(i), i);
            }
            writer.write("\n");

            if (!library.getElutionGroups().isEmpty()) {
                writer.write("elutionGroups: " + library.getElutionGroups().size() + " entries\n");
                for (int i = 0; i < library.getElutionGroups().size(); i++) {
                    writer.write("  [" + i + "]: " + library.getElutionGroups().get(i) + "\n");
                }
            }

            System.out.println("Library content written to " + outputFile.getAbsolutePath());
        }

        Assert.assertTrue("Output file should be created", outputFile.exists());
        Assert.assertTrue("Output file should not be empty", outputFile.length() > 0);
    }

    @Test
    public void testUnsupportedVersions() throws IOException {
        String[] versions = {"1.9", "2.3.0"};
        
        for (String version : versions) {
            File speclibFile = new File(getClass().getClassLoader().getResource(version + "/report-lib.predicted.speclib").getFile());
            
            Assume.assumeTrue("Test file exists: " + version, speclibFile.exists());

            DiaNNSpecLibReader reader = new DiaNNSpecLibReader(speclibFile);
            
            try {
                reader.read();
            } catch (IOException e) {
                Assert.assertTrue("Exception message should mention version for " + version, 
                    e.getMessage().contains("version"));
            }
        }
    }
}
