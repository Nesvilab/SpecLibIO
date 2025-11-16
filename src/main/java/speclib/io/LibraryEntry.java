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

public class LibraryEntry {
    private Peptide target;
    private Peptide decoy;
    private int entryFlags = 0;
    private int proteotypic = 0;
    private String name = "";  // precursor id
    private int pidIndex = 0;
    private int pgIndex = 0;
    private int egId = 0;
    private int bestRun = -1;
    private int peak = 0;
    private int apex = 0;
    private int window = 0;
    private float qvalue = 0.0f;
    private float pgQvalue = 0.0f;
    private float bestFrMz = 0.0f;
    private float ptmQvalue = 0.0f;
    private float siteConf = 0.0f;

    public LibraryEntry() {
        this.target = new Peptide();
    }

    public Peptide getTarget() {
        return target;
    }

    public Peptide getDecoy() {
        return decoy;
    }

    public int getEntryFlags() {
        return entryFlags;
    }

    public int getProteotypic() {
        return proteotypic;
    }

    public String getName() {
        return name;
    }

    public int getPidIndex() {
        return pidIndex;
    }

    public int getPgIndex() {
        return pgIndex;
    }

    public int getEgId() {
        return egId;
    }

    public int getBestRun() {
        return bestRun;
    }

    public int getPeak() {
        return peak;
    }

    public int getApex() {
        return apex;
    }

    public int getWindow() {
        return window;
    }

    public float getQvalue() {
        return qvalue;
    }

    public float getPgQvalue() {
        return pgQvalue;
    }

    public float getBestFrMz() {
        return bestFrMz;
    }

    public float getPtmQvalue() {
        return ptmQvalue;
    }

    public float getSiteConf() {
        return siteConf;
    }

    public void setTarget(Peptide target) {
        this.target = target;
    }

    public void setDecoy(Peptide decoy) {
        this.decoy = decoy;
    }

    public void setEntryFlags(int entryFlags) {
        this.entryFlags = entryFlags;
    }

    public void setProteotypic(int proteotypic) {
        this.proteotypic = proteotypic;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setPidIndex(int pidIndex) {
        this.pidIndex = pidIndex;
    }

    public void setPgIndex(int pgIndex) {
        this.pgIndex = pgIndex;
    }

    public void setEgId(int egId) {
        this.egId = egId;
    }

    public void setBestRun(int bestRun) {
        this.bestRun = bestRun;
    }

    public void setPeak(int peak) {
        this.peak = peak;
    }

    public void setApex(int apex) {
        this.apex = apex;
    }

    public void setWindow(int window) {
        this.window = window;
    }

    public void setQvalue(float qvalue) {
        this.qvalue = qvalue;
    }

    public void setPgQvalue(float pgQvalue) {
        this.pgQvalue = pgQvalue;
    }

    public void setBestFrMz(float bestFrMz) {
        this.bestFrMz = bestFrMz;
    }

    public void setPtmQvalue(float ptmQvalue) {
        this.ptmQvalue = ptmQvalue;
    }

    public void setSiteConf(float siteConf) {
        this.siteConf = siteConf;
    }

    public void write(OutputStream out, int version) throws IOException {
        target.write(out, version);

        ByteBuffer buffer = ByteBuffer.allocate(4);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(decoy != null ? 1 : 0);
        out.write(buffer.array());

        if (decoy != null) {
            decoy.write(out, version);
        }

        buffer.clear();
        buffer.putInt(entryFlags);
        out.write(buffer.array());

        buffer.clear();
        buffer.putInt(proteotypic);
        out.write(buffer.array());

        buffer.clear();
        buffer.putInt(pidIndex);
        out.write(buffer.array());

        writeString(out, name);

        if (version <= -3) {
            buffer.clear();
            buffer.putFloat(pgQvalue);
            out.write(buffer.array());

            buffer.clear();
            buffer.putFloat(ptmQvalue);
            out.write(buffer.array());

            buffer.clear();
            buffer.putFloat(siteConf);
            out.write(buffer.array());
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

    public static LibraryEntry read(InputStream in, int version) throws IOException {
        LibraryEntry entry = new LibraryEntry();
        byte[] buffer = new byte[16];

        entry.target = Peptide.read(in, version);

        int bytesRead = in.read(buffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading LibraryEntry decoy flag");
        }
        int dc = readIntFromBytes(buffer, 0);

        if (dc != 0) {
            entry.decoy = Peptide.read(in, version);
        }

        bytesRead = in.read(buffer, 0, 12);
        if (bytesRead != 12) {
            throw new IOException("Unexpected end of stream while reading LibraryEntry flags");
        }
        entry.entryFlags = readIntFromBytes(buffer, 0);
        entry.proteotypic = readIntFromBytes(buffer, 4);
        entry.pidIndex = readIntFromBytes(buffer, 8);

        entry.name = readString(in);

        if (version <= -3) {
            bytesRead = in.read(buffer, 0, 12);
            if (bytesRead != 12) {
                throw new IOException("Unexpected end of stream while reading LibraryEntry version -3 fields");
            }
            entry.pgQvalue = readFloatFromBytes(buffer, 0);
            entry.ptmQvalue = readFloatFromBytes(buffer, 4);
            entry.siteConf = readFloatFromBytes(buffer, 8);
        }

        return entry;
    }

    private static int readIntFromBytes(byte[] buffer, int offset) {
        return ((buffer[offset] & 0xFF)) |
               ((buffer[offset + 1] & 0xFF) << 8) |
               ((buffer[offset + 2] & 0xFF) << 16) |
               ((buffer[offset + 3] & 0xFF) << 24);
    }

    private static float readFloatFromBytes(byte[] buffer, int offset) {
        int bits = readIntFromBytes(buffer, offset);
        return Float.intBitsToFloat(bits);
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

