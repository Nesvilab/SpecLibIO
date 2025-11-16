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

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class DiaNNSpecLibWriter {
    public static final int LATEST_SUPPORTED_VERSION = -3;

    private final SpectralLibrary library;

    public DiaNNSpecLibWriter(SpectralLibrary library) {
        if (library == null) {
            throw new IllegalArgumentException("Library cannot be null");
        }
        this.library = library;
    }

    public void write(String filePath) throws IOException {
        write(filePath, LATEST_SUPPORTED_VERSION);
    }

    public void write(String filePath, int version) throws IOException {
        if (filePath == null || filePath.isEmpty()) {
            throw new IllegalArgumentException("File path cannot be null or empty");
        }
        if (version > 0) {
            throw new IllegalArgumentException("Version must be negative (e.g., -3)");
        }

        try (FileOutputStream fos = new FileOutputStream(filePath);
             BufferedOutputStream bos = new BufferedOutputStream(fos, 65536)) {
            write(bos, version);
            bos.flush();
        }
    }

    public void write(OutputStream out) throws IOException {
        write(out, LATEST_SUPPORTED_VERSION);
    }

    public void write(OutputStream out, int version) throws IOException {
        if (out == null) {
            throw new IllegalArgumentException("Output stream cannot be null");
        }
        if (version > 0) {
            throw new IllegalArgumentException("Version must be negative (e.g., -3)");
        }

        validateLibrary();
        library.write(out, version);
    }

    public SpectralLibrary getLibrary() {
        return library;
    }

    private void validateLibrary() {
        if (library.getEntries().size() != library.getPrecursors().size()) {
            throw new IllegalStateException(
                String.format("Entries count (%d) must match precursors count (%d)",
                    library.getEntries().size(), library.getPrecursors().size()));
        }

        for (int i = 0; i < library.getEntries().size(); i++) {
            LibraryEntry entry = library.getEntries().get(i);
            if (entry.getTarget() == null) {
                throw new IllegalStateException(
                    String.format("Entry at index %d has null target peptide", i));
            }
            
            String precursor = library.getPrecursors().get(i);
            if (!entry.getName().equals(precursor)) {
                throw new IllegalStateException(
                    String.format("Entry name '%s' at index %d does not match precursor '%s'",
                        entry.getName(), i, precursor));
            }
        }

        if (!library.getElutionGroups().isEmpty() && 
            library.getElutionGroups().size() != library.getEntries().size()) {
            throw new IllegalStateException(
                String.format("Elution groups count (%d) must match entries count (%d) when present",
                    library.getElutionGroups().size(), library.getEntries().size()));
        }
    }

    public static SpectralLibrary createMinimalLibrary(String libraryName) {
        SpectralLibrary library = new SpectralLibrary();
        library.setName(libraryName);
        library.setGenDecoys(false);
        library.setGenCharges(false);
        library.setInferProteotypicity(false);
        return library;
    }

    public static SpectralLibrary createSampleLibrary(String libraryName, 
                                                      String peptideSequence, 
                                                      int charge, 
                                                      float mz) {
        SpectralLibrary library = createMinimalLibrary(libraryName);
        
        LibraryEntry entry = new LibraryEntry();
        String precursorId = peptideSequence + "/" + charge;
        entry.setName(precursorId);
        
        Peptide target = new Peptide();
        target.setIndex(0);
        target.setCharge(charge);
        target.setLength(peptideSequence.length());
        target.setMz(mz);
        target.setiRT(30.0f);
        
        Product fragment1 = new Product(200.1f, 1000.0f, 1);
        Product fragment2 = new Product(300.2f, 2000.0f, 1);
        target.getFragments().add(fragment1);
        target.getFragments().add(fragment2);
        
        entry.setTarget(target);
        
        library.getEntries().add(entry);
        library.getPrecursors().add(precursorId);
        
        return library;
    }
}

