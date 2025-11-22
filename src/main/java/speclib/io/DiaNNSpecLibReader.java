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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

public class DiaNNSpecLibReader {
    public static final int LATEST_SUPPORTED_VERSION = -3;

    private final InputStream inputStream;
    private final boolean ownsStream;

    public DiaNNSpecLibReader(File file) throws IOException {
        this(new BufferedInputStream(new FileInputStream(file), 65536), true);
    }

    public DiaNNSpecLibReader(String filePath) throws IOException {
        this(new File(filePath));
    }

    public DiaNNSpecLibReader(InputStream inputStream) {
        this(inputStream, false);
    }

    private DiaNNSpecLibReader(InputStream inputStream, boolean ownsStream) {
        this.inputStream = inputStream;
        this.ownsStream = ownsStream;
    }

    public SpectralLibrary read() throws IOException {
        try {
            return readLibrary();
        } finally {
            if (ownsStream) {
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private SpectralLibrary readLibrary() throws IOException {
        SpectralLibrary library = new SpectralLibrary();
        byte[] buffer = new byte[16];

        int bytesRead = inputStream.read(buffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading version");
        }
        int version = readIntFromBytes(buffer, 0);

        int gd, gc, ip;
        if (version >= 0) {
            gd = version;
            version = 0;
        } else {
            bytesRead = inputStream.read(buffer, 0, 4);
            if (bytesRead != 4) {
                throw new IOException("Unexpected end of stream while reading gd");
            }
            gd = readIntFromBytes(buffer, 0);
        }

        if (version < LATEST_SUPPORTED_VERSION) {
            throw new IOException(String.format(
                "speclib file has version %d, but this reader only supports up to version %d",
                -1 * version, -1 * LATEST_SUPPORTED_VERSION));
        }

        bytesRead = inputStream.read(buffer, 0, 8);
        if (bytesRead != 8) {
            throw new IOException("Unexpected end of stream while reading flags");
        }
        gc = readIntFromBytes(buffer, 0);
        ip = readIntFromBytes(buffer, 4);

        library.setGenDecoys(gd != 0);
        library.setGenCharges(gc != 0);
        library.setInferProteotypicity(ip != 0);

        library.setName(readString(inputStream));
        library.setFastaNames(readString(inputStream));

        library.setProteins(readIsoformArray(inputStream, version));
        library.setProteinIds(readProteinGroupArray(inputStream, version));

        library.setPrecursors(readStringArray(inputStream));
        library.setNames(readStringArray(inputStream));
        library.setGenes(readStringArray(inputStream));

        bytesRead = inputStream.read(buffer, 0, 16);
        if (bytesRead != 16) {
            throw new IOException("Unexpected end of stream while reading iRT range");
        }
        library.setiRTMin(readDoubleFromBytes(buffer, 0));
        library.setiRTMax(readDoubleFromBytes(buffer, 8));

        library.setEntries(readLibraryEntryArray(inputStream, version));

        if (library.getPrecursors().size() != library.getEntries().size()) {
            throw new IOException(String.format(
                "Precursor count (%d) does not match entry count (%d)",
                library.getPrecursors().size(), library.getEntries().size()));
        }

        for (int i = 0; i < library.getEntries().size(); i++) {
            String entryName = library.getEntries().get(i).getName();
            String precursor = library.getPrecursors().get(i);
            if (!entryName.equals(precursor)) {
                throw new IOException(String.format(
                    "Precursor mismatch between %s and %s at index %d",
                    entryName, precursor, i));
            }
        }

        if (version <= -1) {
            try {
                int peek = inputStream.read();
                if (peek != -1) {
                    byte[] sizeBytes = new byte[4];
                    sizeBytes[0] = (byte) peek;
                    bytesRead = inputStream.read(sizeBytes, 1, 3);
                    if (bytesRead == 3) {
                        int size = readIntFromBytes(sizeBytes, 0);
                        if (size > 0 && size == library.getEntries().size()) {
                            library.setElutionGroups(readIntVector(inputStream, size));
                        }
                    }
                }
            } catch (IOException e) {
                throw new IOException("Failed to read optional elution groups data: " + e.getMessage(), e);
            }
        }

        return library;
    }

    private static int readIntFromBytes(byte[] buffer, int offset) {
        return ((buffer[offset] & 0xFF)) |
               ((buffer[offset + 1] & 0xFF) << 8) |
               ((buffer[offset + 2] & 0xFF) << 16) |
               ((buffer[offset + 3] & 0xFF) << 24);
    }

    private static long readLongFromBytes(byte[] buffer, int offset) {
        return ((buffer[offset] & 0xFFL)) |
               ((buffer[offset + 1] & 0xFFL) << 8) |
               ((buffer[offset + 2] & 0xFFL) << 16) |
               ((buffer[offset + 3] & 0xFFL) << 24) |
               ((buffer[offset + 4] & 0xFFL) << 32) |
               ((buffer[offset + 5] & 0xFFL) << 40) |
               ((buffer[offset + 6] & 0xFFL) << 48) |
               ((buffer[offset + 7] & 0xFFL) << 56);
    }

    private static double readDoubleFromBytes(byte[] buffer, int offset) {
        return Double.longBitsToDouble(readLongFromBytes(buffer, offset));
    }

    private static String readString(InputStream in) throws IOException {
        byte[] sizeBuffer = new byte[4];
        int bytesRead = in.read(sizeBuffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading string length");
        }
        int length = readIntFromBytes(sizeBuffer, 0);

        if (length == 0) {
            return "";
        }

        byte[] strBytes = new byte[length];
        bytesRead = in.read(strBytes, 0, length);
        if (bytesRead != length) {
            throw new IOException("Unexpected end of stream while reading string data");
        }
        return new String(strBytes, StandardCharsets.UTF_8);
    }

    private static List<String> readStringArray(InputStream in) throws IOException {
        byte[] sizeBuffer = new byte[4];
        int bytesRead = in.read(sizeBuffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading string array size");
        }
        int size = readIntFromBytes(sizeBuffer, 0);

        List<String> strings = new ArrayList<>(size);
        for (int i = 0; i < size; i++) {
            strings.add(readString(in));
        }

        return strings;
    }

    private static List<Integer> readIntVector(InputStream in, int size) throws IOException {
        List<Integer> vec = new ArrayList<>(size);
        if (size > 0) {
            byte[] dataBuffer = new byte[size * 4];
            int bytesRead = in.read(dataBuffer, 0, size * 4);
            if (bytesRead != size * 4) {
                throw new IOException("Unexpected end of stream while reading int vector data");
            }

            for (int i = 0; i < size; i++) {
                vec.add(readIntFromBytes(dataBuffer, i * 4));
            }
        }

        return vec;
    }

    private static List<Isoform> readIsoformArray(InputStream in, int version) throws IOException {
        byte[] sizeBuffer = new byte[4];
        int bytesRead = in.read(sizeBuffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading Isoform array size");
        }
        int size = readIntFromBytes(sizeBuffer, 0);

        List<Isoform> isoforms = new ArrayList<>(size);
        for (int i = 0; i < size; i++) {
            isoforms.add(Isoform.read(in, version));
        }

        return isoforms;
    }

    private static List<ProteinGroup> readProteinGroupArray(InputStream in, int version) throws IOException {
        byte[] sizeBuffer = new byte[4];
        int bytesRead = in.read(sizeBuffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading ProteinGroup array size");
        }
        int size = readIntFromBytes(sizeBuffer, 0);

        List<ProteinGroup> proteinGroups = new ArrayList<>(size);
        for (int i = 0; i < size; i++) {
            proteinGroups.add(ProteinGroup.read(in, version));
        }

        return proteinGroups;
    }

    private static List<LibraryEntry> readLibraryEntryArray(InputStream in, int version) throws IOException {
        byte[] sizeBuffer = new byte[4];
        int bytesRead = in.read(sizeBuffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading LibraryEntry array size");
        }
        int size = readIntFromBytes(sizeBuffer, 0);

        List<LibraryEntry> entries = new ArrayList<>(size);
        for (int i = 0; i < size; i++) {
            entries.add(LibraryEntry.read(in, version));
        }

        return entries;
    }
}

