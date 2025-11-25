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
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class DiaNNSpecLibWriterTest {

    @Rule
    public TemporaryFolder tempFolder = new TemporaryFolder();

    @Test
    public void testProductIonCode() {
        Product product = new Product(200.0f, 1000.0f, 1, 2, 5, 0);
        int ionCode = product.getIonCode();
        Assert.assertEquals(24606, ionCode);
    }

    @Test
    public void testBinaryFormatVersion() throws IOException {
        SpectralLibrary library = DiaNNSpecLibWriter.createMinimalLibrary("VersionTest");
        DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
        
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        writer.write(baos, -3);
        
        byte[] data = baos.toByteArray();
        ByteBuffer buffer = ByteBuffer.wrap(data);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        
        Assert.assertEquals(-3, buffer.getInt());
        Assert.assertEquals(0, buffer.getInt());
        Assert.assertEquals(0, buffer.getInt());
        Assert.assertEquals(0, buffer.getInt());
    }

    @Test
    public void testCompleteLibraryStructure() throws IOException {
        SpectralLibrary library = new SpectralLibrary();
        library.setName("CompleteTestLibrary");
        library.setFastaNames("uniprot_human.fasta");
        library.setGenDecoys(false);
        library.setGenCharges(false);
        library.setInferProteotypicity(true);
        library.setiRTMin(5.0);
        library.setiRTMax(95.0);

        Isoform isoform = new Isoform("P12345", "TEST_PROTEIN", "TEST_GENE", "Test protein description", true);
        isoform.setNameIndex(0);
        isoform.setGeneIndex(0);
        isoform.getPrecursors().add(0);
        isoform.getPrecursors().add(1);
        library.getIsoforms().add(isoform);

        ProteinGroup proteinGroup = new ProteinGroup("P12345");
        proteinGroup.setNames("TEST_PROTEIN");
        proteinGroup.setGenes("TEST_GENE");
        proteinGroup.getPrecursors().add(0);
        proteinGroup.getPrecursors().add(1);
        proteinGroup.getIsoforms().add(0);
        library.getProteinGroups().add(proteinGroup);

        library.getNames().add("TEST_PROTEIN");
        library.getGenes().add("TEST_GENE");

        String[] peptides = {"AAAAADLANR", "AAAFEQLQK", "AAANEQLTR"};
        float[] mzs = {472.25f, 503.27f, 487.26f};
        int[] charges = {2, 2, 2};
        float[] rts = {11.12f, 20.11f, 7.85f};

        for (int i = 0; i < peptides.length; i++) {
            LibraryEntry entry = new LibraryEntry();
            String precursorId = peptides[i] + charges[i];
            entry.setName(precursorId);
            entry.setPgIndex(0);

            Precursor target = new Precursor();
            target.setIndex(i);
            target.setCharge(charges[i]);
            target.setLength(peptides[i].length());
            target.setMz(mzs[i]);
            target.setiRT(rts[i]);
            target.setsRT(0.5f);

            for (int j = 0; j < 6; j++) {
                float fragMz = 200.0f + j * 100.0f;
                float fragIntensity = 1000.0f * (6 - j);
                Product fragment = new Product(fragMz, fragIntensity, 1);
                target.getFragments().add(fragment);
            }

            entry.setTarget(target);
            library.getEntries().add(entry);
            library.getPrecursors().add(precursorId);
        }

        DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
        File outputFile = tempFolder.newFile("complete-library.speclib");
        writer.write(outputFile.getAbsolutePath());

        Assert.assertTrue(outputFile.exists());

        try (FileInputStream fis = new FileInputStream(outputFile)) {
            byte[] header = new byte[16];
            fis.read(header);
            ByteBuffer buffer = ByteBuffer.wrap(header);
            buffer.order(ByteOrder.LITTLE_ENDIAN);

            Assert.assertEquals("Version should be -3", -3, buffer.getInt());
            Assert.assertEquals("genDecoys should be false", 0, buffer.getInt());
            Assert.assertEquals("genCharges should be false", 0, buffer.getInt());
            Assert.assertEquals("inferProteotypicity should be true", 1, buffer.getInt());
        }
    }

    @Test
    public void testDeterministicOutput() throws IOException {
        SpectralLibrary library = DiaNNSpecLibWriter.createSampleLibrary(
            "DeterministicTest", "PEPTIDEK", 2, 450.25f);

        DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);

        ByteArrayOutputStream baos1 = new ByteArrayOutputStream();
        writer.write(baos1);
        byte[] data1 = baos1.toByteArray();

        ByteArrayOutputStream baos2 = new ByteArrayOutputStream();
        writer.write(baos2);
        byte[] data2 = baos2.toByteArray();

        Assert.assertArrayEquals("Output should be deterministic", data1, data2);
    }

    @Test
    public void testRoundTrip_DiannSwath() throws IOException {
        File originalFile = new File(getClass().getClassLoader().getResource("diann-swath.speclib").getFile());
        org.junit.Assume.assumeTrue("Test file exists", originalFile.exists());
        
        testRoundTripForFile(originalFile, "diann-swath-roundtrip.speclib");
    }

    @Test
    public void testRoundTrip_DiannModTest() throws IOException {
        File originalFile = new File(getClass().getClassLoader().getResource("diann-mod-test.tsv.speclib").getFile());
        org.junit.Assume.assumeTrue("Test file exists", originalFile.exists());
        
        testRoundTripForFile(originalFile, "diann-mod-test-roundtrip.speclib");
    }

    @Test
    public void testRoundTrip_DiannMassMods() throws IOException {
        File originalFile = new File(getClass().getClassLoader().getResource("diann-mass-mods.tsv.speclib").getFile());
        org.junit.Assume.assumeTrue("Test file exists", originalFile.exists());
        
        testRoundTripForFile(originalFile, "diann-mass-mods-roundtrip.speclib");
    }

    @Test
    public void testRoundTrip_DiannDiaPasef() throws IOException {
        File originalFile = new File(getClass().getClassLoader().getResource("diann-hela-diapasef-lib.speclib").getFile());
        org.junit.Assume.assumeTrue("Test file exists", originalFile.exists());
        
        testRoundTripForFile(originalFile, "diann-hela-diapasef-roundtrip.speclib");
    }

    @Test
    public void testRoundTrip_LibraryTsv() throws IOException {
        File originalFile = new File(getClass().getClassLoader().getResource("library.tsv.speclib").getFile());
        org.junit.Assume.assumeTrue("Test file exists", originalFile.exists());
        
        testRoundTripForFile(originalFile, "library-tsv-roundtrip.speclib");
    }

    private void testRoundTripForFile(File originalFile, String newFileName) throws IOException {
        DiaNNSpecLibReader reader = new DiaNNSpecLibReader(originalFile);
        SpectralLibrary originalLibrary = reader.read();

        File outputFile = tempFolder.newFile(newFileName);
        DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(originalLibrary);
        writer.write(outputFile.getAbsolutePath());

        Assert.assertTrue("Written file should exist", outputFile.exists());
        Assert.assertTrue("Written file should have data", outputFile.length() > 0);

        byte[] originalBytes = readFileBytes(originalFile);
        byte[] writtenBytes = readFileBytes(outputFile);
        Assert.assertEquals("File sizes should match", originalBytes.length, writtenBytes.length);
        Assert.assertArrayEquals("Binary content should match", originalBytes, writtenBytes);

        DiaNNSpecLibReader reReader = new DiaNNSpecLibReader(outputFile);
        SpectralLibrary rereadLibrary = reReader.read();

        compareLibraries(originalLibrary, rereadLibrary, originalFile.getName());
    }

    private byte[] readFileBytes(File file) throws IOException {
        try (FileInputStream fis = new FileInputStream(file)) {
            byte[] buffer = new byte[(int) file.length()];
            int bytesRead = fis.read(buffer);
            if (bytesRead != buffer.length) {
                throw new IOException("Could not read entire file");
            }
            return buffer;
        }
    }

    private void compareLibraries(SpectralLibrary original, SpectralLibrary reread, String fileName) {
        Assert.assertEquals("Library name mismatch for " + fileName, 
            original.getName(), reread.getName());
        Assert.assertEquals("FASTA names mismatch for " + fileName, 
            original.getFastaNames(), reread.getFastaNames());
        Assert.assertEquals("GenDecoys flag mismatch for " + fileName, 
            original.isGenDecoys(), reread.isGenDecoys());
        Assert.assertEquals("GenCharges flag mismatch for " + fileName, 
            original.isGenCharges(), reread.isGenCharges());
        Assert.assertEquals("InferProteotypicity flag mismatch for " + fileName, 
            original.isInferProteotypicity(), reread.isInferProteotypicity());
        Assert.assertEquals("iRT min mismatch for " + fileName, 
            original.getiRTMin(), reread.getiRTMin(), 0.0001);
        Assert.assertEquals("iRT max mismatch for " + fileName, 
            original.getiRTMax(), reread.getiRTMax(), 0.0001);

        Assert.assertEquals("Precursors size mismatch for " + fileName, 
            original.getPrecursors().size(), reread.getPrecursors().size());
        for (int i = 0; i < original.getPrecursors().size(); i++) {
            Assert.assertEquals("Precursor[" + i + "] mismatch for " + fileName, 
                original.getPrecursors().get(i), reread.getPrecursors().get(i));
        }

        Assert.assertEquals("Names size mismatch for " + fileName, 
            original.getNames().size(), reread.getNames().size());
        for (int i = 0; i < original.getNames().size(); i++) {
            Assert.assertEquals("Name[" + i + "] mismatch for " + fileName, 
                original.getNames().get(i), reread.getNames().get(i));
        }

        Assert.assertEquals("Genes size mismatch for " + fileName, 
            original.getGenes().size(), reread.getGenes().size());
        for (int i = 0; i < original.getGenes().size(); i++) {
            Assert.assertEquals("Gene[" + i + "] mismatch for " + fileName, 
                original.getGenes().get(i), reread.getGenes().get(i));
        }

        Assert.assertEquals("Isoforms size mismatch for " + fileName, 
            original.getIsoforms().size(), reread.getIsoforms().size());
        for (int i = 0; i < original.getIsoforms().size(); i++) {
            compareIsoforms(original.getIsoforms().get(i), reread.getIsoforms().get(i), fileName, i);
        }

        Assert.assertEquals("ProteinGroups size mismatch for " + fileName, 
            original.getProteinGroups().size(), reread.getProteinGroups().size());
        for (int i = 0; i < original.getProteinGroups().size(); i++) {
            compareProteinGroups(original.getProteinGroups().get(i), reread.getProteinGroups().get(i), fileName, i);
        }

        Assert.assertEquals("Entries size mismatch for " + fileName, 
            original.getEntries().size(), reread.getEntries().size());
        for (int i = 0; i < original.getEntries().size(); i++) {
            compareLibraryEntries(original.getEntries().get(i), reread.getEntries().get(i), fileName, i);
        }

        Assert.assertEquals("ElutionGroups size mismatch for " + fileName, 
            original.getElutionGroups().size(), reread.getElutionGroups().size());
        for (int i = 0; i < original.getElutionGroups().size(); i++) {
            Assert.assertEquals("ElutionGroup[" + i + "] mismatch for " + fileName, 
                original.getElutionGroups().get(i), reread.getElutionGroups().get(i));
        }
    }

    private void compareIsoforms(Isoform original, Isoform reread, String fileName, int index) {
        Assert.assertEquals("Isoform[" + index + "] id mismatch for " + fileName, 
            original.getId(), reread.getId());
        Assert.assertEquals("Isoform[" + index + "] name mismatch for " + fileName, 
            original.getName(), reread.getName());
        Assert.assertEquals("Isoform[" + index + "] gene mismatch for " + fileName, 
            original.getGene(), reread.getGene());
        Assert.assertEquals("Isoform[" + index + "] description mismatch for " + fileName, 
            original.getDescription(), reread.getDescription());
        Assert.assertEquals("Isoform[" + index + "] swissprot flag mismatch for " + fileName, 
            original.isSwissprot(), reread.isSwissprot());
        Assert.assertEquals("Isoform[" + index + "] nameIndex mismatch for " + fileName, 
            original.getNameIndex(), reread.getNameIndex());
        Assert.assertEquals("Isoform[" + index + "] geneIndex mismatch for " + fileName, 
            original.getGeneIndex(), reread.getGeneIndex());
        Assert.assertEquals("Isoform[" + index + "] precursors size mismatch for " + fileName, 
            original.getPrecursors().size(), reread.getPrecursors().size());
        Assert.assertEquals("Isoform[" + index + "] precursors mismatch for " + fileName, 
            original.getPrecursors(), reread.getPrecursors());
    }

    private void compareProteinGroups(ProteinGroup original, ProteinGroup reread, String fileName, int index) {
        Assert.assertEquals("ProteinGroup[" + index + "] ids mismatch for " + fileName, 
            original.getIds(), reread.getIds());
        Assert.assertEquals("ProteinGroup[" + index + "] names mismatch for " + fileName, 
            original.getNames(), reread.getNames());
        Assert.assertEquals("ProteinGroup[" + index + "] genes mismatch for " + fileName, 
            original.getGenes(), reread.getGenes());
        Assert.assertEquals("ProteinGroup[" + index + "] precursors size mismatch for " + fileName, 
            original.getPrecursors().size(), reread.getPrecursors().size());
        for (int i = 0; i < original.getPrecursors().size(); i++) {
            Assert.assertEquals("ProteinGroup[" + index + "] precursor[" + i + "] mismatch for " + fileName, 
                original.getPrecursors().get(i), reread.getPrecursors().get(i));
        }
        Assert.assertEquals("ProteinGroup[" + index + "] isoforms size mismatch for " + fileName, 
            original.getIsoforms().size(), reread.getIsoforms().size());
        Assert.assertEquals("ProteinGroup[" + index + "] isoforms mismatch for " + fileName, 
            original.getIsoforms(), reread.getIsoforms());
        Assert.assertEquals("ProteinGroup[" + index + "] nameIndices size mismatch for " + fileName, 
            original.getNameIndices().size(), reread.getNameIndices().size());
        for (int i = 0; i < original.getNameIndices().size(); i++) {
            Assert.assertEquals("ProteinGroup[" + index + "] nameIndex[" + i + "] mismatch for " + fileName, 
                original.getNameIndices().get(i), reread.getNameIndices().get(i));
        }
        Assert.assertEquals("ProteinGroup[" + index + "] geneIndices size mismatch for " + fileName, 
            original.getGeneIndices().size(), reread.getGeneIndices().size());
        for (int i = 0; i < original.getGeneIndices().size(); i++) {
            Assert.assertEquals("ProteinGroup[" + index + "] geneIndex[" + i + "] mismatch for " + fileName, 
                original.getGeneIndices().get(i), reread.getGeneIndices().get(i));
        }
    }

    private void compareLibraryEntries(LibraryEntry original, LibraryEntry reread, String fileName, int index) {
        Assert.assertEquals("Entry[" + index + "] name mismatch for " + fileName, 
            original.getName(), reread.getName());
        Assert.assertEquals("Entry[" + index + "] pidIndex mismatch for " + fileName, 
            original.getPgIndex(), reread.getPgIndex());
        Assert.assertEquals("Entry[" + index + "] proteotypic mismatch for " + fileName, 
            original.getProteotypic(), reread.getProteotypic());
        Assert.assertEquals("Entry[" + index + "] pgQvalue mismatch for " + fileName, 
            original.getPgQvalue(), reread.getPgQvalue(), 0.0001f);
        Assert.assertEquals("Entry[" + index + "] ptmQvalue mismatch for " + fileName, 
            original.getPtmQvalue(), reread.getPtmQvalue(), 0.0001f);
        Assert.assertEquals("Entry[" + index + "] siteConf mismatch for " + fileName, 
            original.getSiteConf(), reread.getSiteConf(), 0.0001f);

        comparePrecursors(original.getTarget(), reread.getTarget(), fileName, index, "target");
        
        if (original.getDecoy() != null || reread.getDecoy() != null) {
            Assert.assertNotNull("Entry[" + index + "] original decoy should not be null for " + fileName, 
                original.getDecoy());
            Assert.assertNotNull("Entry[" + index + "] reread decoy should not be null for " + fileName, 
                reread.getDecoy());
            comparePrecursors(original.getDecoy(), reread.getDecoy(), fileName, index, "decoy");
        }
    }

    private void comparePrecursors(Precursor original, Precursor reread, String fileName, int entryIndex, String type) {
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " index mismatch for " + fileName, 
            original.getIndex(), reread.getIndex());
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " charge mismatch for " + fileName, 
            original.getCharge(), reread.getCharge());
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " mz mismatch for " + fileName, 
            original.getMz(), reread.getMz(), 0.0001f);
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " length mismatch for " + fileName, 
            original.getLength(), reread.getLength());
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " iRT mismatch for " + fileName, 
            original.getiRT(), reread.getiRT(), 0.001f);
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " sRT mismatch for " + fileName, 
            original.getsRT(), reread.getsRT(), 0.001f);
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " iIM mismatch for " + fileName, 
            original.getiIM(), reread.getiIM(), 0.0001f);
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " sIM mismatch for " + fileName, 
            original.getsIM(), reread.getsIM(), 0.0001f);
        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " libQvalue mismatch for " + fileName, 
            original.getLibQvalue(), reread.getLibQvalue(), 0.0001f);

        Assert.assertEquals("Entry[" + entryIndex + "]." + type + " fragments size mismatch for " + fileName, 
            original.getFragments().size(), reread.getFragments().size());
        for (int i = 0; i < original.getFragments().size(); i++) {
            compareProducts(original.getFragments().get(i), reread.getFragments().get(i), 
                fileName, entryIndex, type, i);
        }
    }

    private void compareProducts(Product original, Product reread, String fileName, int entryIndex, 
                                 String peptideType, int fragIndex) {
        Assert.assertEquals("Entry[" + entryIndex + "]." + peptideType + ".fragment[" + fragIndex + "] mz mismatch for " + fileName, 
            original.getMz(), reread.getMz(), 0.0001f);
        Assert.assertEquals("Entry[" + entryIndex + "]." + peptideType + ".fragment[" + fragIndex + "] height mismatch for " + fileName, 
            original.getHeight(), reread.getHeight(), 0.01f);
        Assert.assertEquals("Entry[" + entryIndex + "]." + peptideType + ".fragment[" + fragIndex + "] charge mismatch for " + fileName, 
            original.getCharge(), reread.getCharge());
        Assert.assertEquals("Entry[" + entryIndex + "]." + peptideType + ".fragment[" + fragIndex + "] type mismatch for " + fileName, 
            original.getType(), reread.getType());
        Assert.assertEquals("Entry[" + entryIndex + "]." + peptideType + ".fragment[" + fragIndex + "] index mismatch for " + fileName, 
            original.getIndex(), reread.getIndex());
        Assert.assertEquals("Entry[" + entryIndex + "]." + peptideType + ".fragment[" + fragIndex + "] loss mismatch for " + fileName, 
            original.getLoss(), reread.getLoss());
    }
}

