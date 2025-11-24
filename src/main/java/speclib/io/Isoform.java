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
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashSet;
import java.util.Set;

public class Isoform {
    private String id = "";
    private String name = "";
    private String gene = "";
    private String description = "";
    private Set<Integer> precursors = new LinkedHashSet<>();
    private int nameIndex = 0;
    private int geneIndex = 0;
    private boolean swissprot = false;

    public Isoform() {
    }

    public Isoform(String id) {
        this.id = id;
    }

    public Isoform(String id, String name, String gene, String description, boolean swissprot) {
        this.id = id;
        this.name = name;
        this.gene = gene;
        this.description = description;
        this.swissprot = swissprot;
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public String getGene() {
        return gene;
    }

    public String getDescription() {
        return description;
    }

    public Set<Integer> getPrecursors() {
        return precursors;
    }

    public int getNameIndex() {
        return nameIndex;
    }

    public int getGeneIndex() {
        return geneIndex;
    }

    public boolean isSwissprot() {
        return swissprot;
    }

    public void setId(String id) {
        this.id = id;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setGene(String gene) {
        this.gene = gene;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public void setPrecursors(Set<Integer> precursors) {
        this.precursors = precursors;
    }

    public void setNameIndex(int nameIndex) {
        this.nameIndex = nameIndex;
    }

    public void setGeneIndex(int geneIndex) {
        this.geneIndex = geneIndex;
    }

    public void setSwissprot(boolean swissprot) {
        this.swissprot = swissprot;
    }

    public void write(OutputStream out, int version) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(8);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        buffer.putInt(swissprot ? 1 : 0);
        buffer.putInt(precursors.size());
        out.write(buffer.array(), 0, 8);

        writeString(out, id);
        writeString(out, name);
        writeString(out, gene);

        buffer.clear();
        buffer.putInt(nameIndex);
        buffer.putInt(geneIndex);
        out.write(buffer.array(), 0, 8);

        if (!precursors.isEmpty()) {
            ByteBuffer precursorBuffer = ByteBuffer.allocate(precursors.size() * 4);
            precursorBuffer.order(ByteOrder.LITTLE_ENDIAN);
            for (Integer precursor : precursors) {
                precursorBuffer.putInt(precursor);
            }
            out.write(precursorBuffer.array());
        }
    }

    private void writeString(OutputStream out, String str) throws IOException {
        byte[] bytes = str.getBytes(StandardCharsets.UTF_8);
        ByteBuffer buffer = ByteBuffer.allocate(4);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(bytes.length);
        out.write(buffer.array());
        if (bytes.length > 0) {
            out.write(bytes);
        }
    }

    public static Isoform read(InputStream in, int version) throws IOException {
        Isoform isoform = new Isoform();
        byte[] buffer = new byte[8];

        int bytesRead = in.read(buffer, 0, 8);
        if (bytesRead != 8) {
            throw new IOException("Unexpected end of stream while reading Isoform header");
        }
        int sp = readIntFromBytes(buffer, 0);
        int size = readIntFromBytes(buffer, 4);
        isoform.swissprot = (sp == 1);

        isoform.id = readString(in);
        isoform.name = readString(in);
        isoform.gene = readString(in);

        bytesRead = in.read(buffer, 0, 8);
        if (bytesRead != 8) {
            throw new IOException("Unexpected end of stream while reading Isoform indices");
        }
        isoform.nameIndex = readIntFromBytes(buffer, 0);
        isoform.geneIndex = readIntFromBytes(buffer, 4);

        isoform.precursors.clear();
        for (int i = 0; i < size; i++) {
            bytesRead = in.read(buffer, 0, 4);
            if (bytesRead != 4) {
                throw new IOException("Unexpected end of stream while reading Isoform precursors");
            }
            int pr = readIntFromBytes(buffer, 0);
            if (pr >= 0) {
                isoform.precursors.add(pr);
            }
        }

        return isoform;
    }

    private static int readIntFromBytes(byte[] buffer, int offset) {
        return ((buffer[offset] & 0xFF)) |
               ((buffer[offset + 1] & 0xFF) << 8) |
               ((buffer[offset + 2] & 0xFF) << 16) |
               ((buffer[offset + 3] & 0xFF) << 24);
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
}

