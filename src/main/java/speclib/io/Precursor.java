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
import java.util.ArrayList;
import java.util.List;

public class Precursor {
    private int index = 0;
    private int charge = 0;
    private int length = 0;
    private int noCal = 0;
    private float mz = 0.0f;
    private float iRT = 0.0f;
    private float sRT = 0.0f;
    private float libQvalue = 0.0f;
    private float iIM = 0.0f;
    private float sIM = 0.0f;
    private List<Product> fragments = new ArrayList<>();

    public Precursor() {
    }

    public int getIndex() {
        return index;
    }

    public int getCharge() {
        return charge;
    }

    public int getLength() {
        return length;
    }

    public int getNoCal() {
        return noCal;
    }

    public float getMz() {
        return mz;
    }

    public float getiRT() {
        return iRT;
    }

    public float getsRT() {
        return sRT;
    }

    public float getLibQvalue() {
        return libQvalue;
    }

    public float getiIM() {
        return iIM;
    }

    public float getsIM() {
        return sIM;
    }

    public List<Product> getFragments() {
        return fragments;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public void setCharge(int charge) {
        this.charge = charge;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public void setNoCal(int noCal) {
        this.noCal = noCal;
    }

    public void setMz(float mz) {
        this.mz = mz;
    }

    public void setiRT(float iRT) {
        this.iRT = iRT;
    }

    public void setsRT(float sRT) {
        this.sRT = sRT;
    }

    public void setLibQvalue(float libQvalue) {
        this.libQvalue = libQvalue;
    }

    public void setiIM(float iIM) {
        this.iIM = iIM;
    }

    public void setsIM(float sIM) {
        this.sIM = sIM;
    }

    public void setFragments(List<Product> fragments) {
        this.fragments = fragments;
    }

    public void write(OutputStream out, int version) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(36);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        buffer.putInt(index);
        buffer.putInt(charge);
        buffer.putInt(length);

        buffer.putFloat(mz);
        buffer.putFloat(iRT);
        buffer.putFloat(sRT);

        if (version <= -2) {
            buffer.putFloat(libQvalue);
            buffer.putFloat(iIM);
            buffer.putFloat(sIM);
        }

        out.write(buffer.array(), 0, buffer.position());

        writeFragmentsVector(out);
    }

    private void writeFragmentsVector(OutputStream out) throws IOException {
        ByteBuffer sizeBuffer = ByteBuffer.allocate(4);
        sizeBuffer.order(ByteOrder.LITTLE_ENDIAN);
        sizeBuffer.putInt(fragments.size());
        out.write(sizeBuffer.array());

        if (!fragments.isEmpty()) {
            int productSize = 12;
            ByteBuffer fragmentsBuffer = ByteBuffer.allocate(fragments.size() * productSize);
            fragmentsBuffer.order(ByteOrder.LITTLE_ENDIAN);

            for (Product fragment : fragments) {
                fragmentsBuffer.putFloat(fragment.getMz());
                fragmentsBuffer.putFloat(fragment.getHeight());
                fragmentsBuffer.put(fragment.getCharge());
                fragmentsBuffer.put(fragment.getType());
                fragmentsBuffer.put(fragment.getIndex());
                fragmentsBuffer.put(fragment.getLoss());
            }

            out.write(fragmentsBuffer.array());
        }
    }

    public static Precursor read(InputStream in, int version) throws IOException {
        Precursor precursor = new Precursor();
        byte[] buffer = new byte[48];

        int bytesRead = in.read(buffer, 0, 12);
        if (bytesRead != 12) {
            throw new IOException("Unexpected end of stream while reading Precursor ints");
        }
        precursor.index = readIntFromBytes(buffer, 0);
        precursor.charge = readIntFromBytes(buffer, 4);
        precursor.length = readIntFromBytes(buffer, 8);

        bytesRead = in.read(buffer, 0, 12);
        if (bytesRead != 12) {
            throw new IOException("Unexpected end of stream while reading Precursor floats");
        }
        precursor.mz = readFloatFromBytes(buffer, 0);
        precursor.iRT = readFloatFromBytes(buffer, 4);
        precursor.sRT = readFloatFromBytes(buffer, 8);

        if (version <= -2) {
            bytesRead = in.read(buffer, 0, 12);
            if (bytesRead != 12) {
                throw new IOException("Unexpected end of stream while reading Precursor version -2 fields");
            }
            precursor.libQvalue = readFloatFromBytes(buffer, 0);
            precursor.iIM = readFloatFromBytes(buffer, 4);
            precursor.sIM = readFloatFromBytes(buffer, 8);
        }

        precursor.fragments = readFragmentsVector(in);

        return precursor;
    }

    private static List<Product> readFragmentsVector(InputStream in) throws IOException {
        byte[] buffer = new byte[4];
        int bytesRead = in.read(buffer, 0, 4);
        if (bytesRead != 4) {
            throw new IOException("Unexpected end of stream while reading fragments vector size");
        }
        int size = readIntFromBytes(buffer, 0);

        List<Product> fragments = new ArrayList<>(size);
        if (size > 0) {
            byte[] productBuffer = new byte[12];
            for (int i = 0; i < size; i++) {
                fragments.add(Product.read(productBuffer, in));
            }
        }

        return fragments;
    }

    private static int readIntFromBytes(byte[] buffer, int offset) {
        return ((buffer[offset] & 0xFF)) |
               ((buffer[offset + 1] & 0xFF) << 8) |
               ((buffer[offset + 2] & 0xFF) << 16) |
               ((buffer[offset + 3] & 0xFF) << 24);
    }

    private static float readFloatFromBytes(byte[] buffer, int offset) {
        return Float.intBitsToFloat(readIntFromBytes(buffer, offset));
    }
}

