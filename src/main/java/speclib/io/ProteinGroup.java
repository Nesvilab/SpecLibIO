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
import java.util.ArrayList;
import java.util.List;

public class ProteinGroup implements Comparable<ProteinGroup> {
    private String ids = "";
    private String names = "";
    private String genes = "";
    private List<Integer> precursors = new ArrayList<>();
    private List<Integer> isoforms = new ArrayList<>();
    private List<Integer> nameIndices = new ArrayList<>();
    private List<Integer> geneIndices = new ArrayList<>();

    public ProteinGroup() {
    }

    public ProteinGroup(String ids) {
        this.ids = ids;
    }

    @Override
    public int compareTo(ProteinGroup other) {
        return this.ids.compareTo(other.ids);
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        ProteinGroup that = (ProteinGroup) obj;
        return ids.equals(that.ids);
    }

    @Override
    public int hashCode() {
        return ids.hashCode();
    }

    public String getIds() {
        return ids;
    }

    public String getNames() {
        return names;
    }

    public String getGenes() {
        return genes;
    }

    public List<Integer> getPrecursors() {
        return precursors;
    }

    public List<Integer> getIsoforms() {
        return isoforms;
    }

    public List<Integer> getNameIndices() {
        return nameIndices;
    }

    public List<Integer> getGeneIndices() {
        return geneIndices;
    }

    public void setIds(String ids) {
        this.ids = ids;
    }

    public void setNames(String names) {
        this.names = names;
    }

    public void setGenes(String genes) {
        this.genes = genes;
    }

    public void setPrecursors(List<Integer> precursors) {
        this.precursors = precursors;
    }

    public void setIsoforms(List<Integer> isoforms) {
        this.isoforms = isoforms;
    }

    public void setNameIndices(List<Integer> nameIndices) {
        this.nameIndices = nameIndices;
    }

    public void setGeneIndices(List<Integer> geneIndices) {
        this.geneIndices = geneIndices;
    }

    public void write(OutputStream out, int version) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(4);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        buffer.putInt(isoforms.size());
        out.write(buffer.array());

        writeString(out, ids);
        writeString(out, names);
        writeString(out, genes);

        writeIntVector(out, nameIndices);

        writeIntVector(out, geneIndices);

        writeIntVector(out, precursors);

        if (!isoforms.isEmpty()) {
            ByteBuffer isoformBuffer = ByteBuffer.allocate(isoforms.size() * 4);
            isoformBuffer.order(ByteOrder.LITTLE_ENDIAN);
            for (Integer isoform : isoforms) {
                isoformBuffer.putInt(isoform);
            }
            out.write(isoformBuffer.array());
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

    private void writeIntVector(OutputStream out, List<Integer> vec) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(4);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(vec.size());
        out.write(buffer.array());

        if (!vec.isEmpty()) {
            ByteBuffer dataBuffer = ByteBuffer.allocate(vec.size() * 4);
            dataBuffer.order(ByteOrder.LITTLE_ENDIAN);
            for (Integer val : vec) {
                dataBuffer.putInt(val);
            }
            out.write(dataBuffer.array());
        }
    }

    public static ProteinGroup read(InputStream in, int version) throws IOException {
        ProteinGroup proteinGroup = new ProteinGroup();
        byte[] buffer = new byte[4];

        int bytesRead = in.read(buffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading ProteinGroup size");
        }
        int sizeP = readIntFromBytes(buffer, 0);

        proteinGroup.ids = readString(in);
        proteinGroup.names = readString(in);
        proteinGroup.genes = readString(in);

        proteinGroup.nameIndices = readIntVector(in);
        proteinGroup.geneIndices = readIntVector(in);
        proteinGroup.precursors = readIntVector(in);

        proteinGroup.isoforms.clear();
        for (int i = 0; i < sizeP; i++) {
            bytesRead = in.read(buffer, 0, 4);
            if (bytesRead != 4) {
                throw new IOException("Unexpected end of stream while reading ProteinGroup proteins");
            }
            int p = readIntFromBytes(buffer, 0);
            if (p >= 0) {
                proteinGroup.isoforms.add(p);
            }
        }

        return proteinGroup;
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

    private static List<Integer> readIntVector(InputStream in) throws IOException {
        byte[] sizeBuffer = new byte[4];
        int bytesRead = in.read(sizeBuffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading int vector size");
        }
        int size = readIntFromBytes(sizeBuffer, 0);

        List<Integer> vec = new ArrayList<>(size);
        if (size > 0) {
            byte[] dataBuffer = new byte[size * 4];
            bytesRead = in.read(dataBuffer, 0, size * 4);
            if (bytesRead != size * 4) {
                throw new IOException("Unexpected end of stream while reading int vector data");
            }

            for (int i = 0; i < size; i++) {
                vec.add(readIntFromBytes(dataBuffer, i * 4));
            }
        }

        return vec;
    }
}

