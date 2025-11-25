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
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

public class SpectralLibrary {
    private String name = "";
    private String fastaNames = "";
    private List<Isoform> isoforms = new ArrayList<>();
    private List<ProteinGroup> proteinGroups = new ArrayList<>();
    private List<String> precursors = new ArrayList<>();
    private List<String> names = new ArrayList<>();
    private List<String> genes = new ArrayList<>();

    private int skipped = 0;
    private double iRTMin = 0.0;
    private double iRTMax = 0.0;
    private boolean genDecoys = true;
    private boolean genCharges = true;
    private boolean inferProteotypicity = true;

    private List<LibraryEntry> entries = new ArrayList<>();
    private List<Integer> elutionGroups = new ArrayList<>();

    public SpectralLibrary() {
    }

    public String getName() {
        return name;
    }

    public String getFastaNames() {
        return fastaNames;
    }

    public List<Isoform> getIsoforms() {
        return isoforms;
    }

    public List<ProteinGroup> getProteinGroups() {
        return proteinGroups;
    }

    public List<String> getPrecursors() {
        return precursors;
    }

    public List<String> getNames() {
        return names;
    }

    public List<String> getGenes() {
        return genes;
    }

    public int getSkipped() {
        return skipped;
    }

    public double getiRTMin() {
        return iRTMin;
    }

    public double getiRTMax() {
        return iRTMax;
    }

    public boolean isGenDecoys() {
        return genDecoys;
    }

    public boolean isGenCharges() {
        return genCharges;
    }

    public boolean isInferProteotypicity() {
        return inferProteotypicity;
    }

    public List<LibraryEntry> getEntries() {
        return entries;
    }

    public List<Integer> getElutionGroups() {
        return elutionGroups;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setFastaNames(String fastaNames) {
        this.fastaNames = fastaNames;
    }

    public void setIsoforms(List<Isoform> isoforms) {
        this.isoforms = isoforms;
    }

    public void setProteinGroups(List<ProteinGroup> proteinGroups) {
        this.proteinGroups = proteinGroups;
    }

    public void setPrecursors(List<String> precursors) {
        this.precursors = precursors;
    }

    public void setNames(List<String> names) {
        this.names = names;
    }

    public void setGenes(List<String> genes) {
        this.genes = genes;
    }

    public void setSkipped(int skipped) {
        this.skipped = skipped;
    }

    public void setiRTMin(double iRTMin) {
        this.iRTMin = iRTMin;
    }

    public void setiRTMax(double iRTMax) {
        this.iRTMax = iRTMax;
    }

    public void setGenDecoys(boolean genDecoys) {
        this.genDecoys = genDecoys;
    }

    public void setGenCharges(boolean genCharges) {
        this.genCharges = genCharges;
    }

    public void setInferProteotypicity(boolean inferProteotypicity) {
        this.inferProteotypicity = inferProteotypicity;
    }

    public void setEntries(List<LibraryEntry> entries) {
        this.entries = entries;
    }

    public void setElutionGroups(List<Integer> elutionGroups) {
        this.elutionGroups = elutionGroups;
    }

    public void write(OutputStream out) throws IOException {
        write(out, -3);
    }

    public void write(OutputStream out, int version) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(8);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        if (version >= 0) {
            buffer.putInt(genDecoys ? 1 : 0);
            out.write(buffer.array(), 0, 4);
        } else {
            buffer.putInt(version);
            out.write(buffer.array(), 0, 4);

            buffer.clear();
            buffer.putInt(genDecoys ? 1 : 0);
            out.write(buffer.array(), 0, 4);
        }

        buffer.clear();
        buffer.putInt(genCharges ? 1 : 0);
        out.write(buffer.array(), 0, 4);

        buffer.clear();
        buffer.putInt(inferProteotypicity ? 1 : 0);
        out.write(buffer.array(), 0, 4);

        writeString(out, name);
        writeString(out, fastaNames);

        writeArray(out, isoforms, version);
        writeArray(out, proteinGroups, version);

        writeStrings(out, precursors);
        writeStrings(out, names);
        writeStrings(out, genes);

        buffer.clear();
        buffer.putDouble(iRTMin);
        out.write(buffer.array(), 0, 8);

        buffer.clear();
        buffer.putDouble(iRTMax);
        out.write(buffer.array(), 0, 8);

        writeArray(out, entries, version);

        if (version <= -1 && !elutionGroups.isEmpty()) {
            writeIntVector(out, elutionGroups);
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

    private void writeStrings(OutputStream out, List<String> strings) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(4);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(strings.size());
        out.write(buffer.array());

        for (String str : strings) {
            writeString(out, str);
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

    private <T> void writeArray(OutputStream out, List<T> array, int version) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(4);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(array.size());
        out.write(buffer.array());

        for (T item : array) {
            if (item instanceof Isoform) {
                ((Isoform) item).write(out, version);
            } else if (item instanceof ProteinGroup) {
                ((ProteinGroup) item).write(out, version);
            } else if (item instanceof LibraryEntry) {
                ((LibraryEntry) item).write(out, version);
            }
        }
    }
}

