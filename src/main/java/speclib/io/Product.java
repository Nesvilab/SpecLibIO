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

public class Product {
    private float mz = 0.0f;
    private float height = 0.0f;
    private byte charge = 0;
    private byte type = 0;
    private byte index = 0;
    private byte loss = 0;

    public Product() {
    }

    public Product(float mz, float height, int charge) {
        this.mz = mz;
        this.height = height;
        this.charge = (byte) charge;
    }

    public Product(float mz, float height, int charge, int type, int index, int loss) {
        this.mz = mz;
        this.height = height;
        this.charge = (byte) charge;
        this.type = (byte) type;
        this.index = (byte) index;
        this.loss = (byte) loss;
    }

    public float getMz() {
        return mz;
    }

    public float getHeight() {
        return height;
    }

    public byte getCharge() {
        return charge;
    }

    public byte getType() {
        return type;
    }

    public byte getIndex() {
        return index;
    }

    public byte getLoss() {
        return loss;
    }

    public void setMz(float mz) {
        this.mz = mz;
    }

    public void setHeight(float height) {
        this.height = height;
    }

    public void setCharge(byte charge) {
        this.charge = charge;
    }

    public void setType(byte type) {
        this.type = type;
    }

    public void setIndex(byte index) {
        this.index = index;
    }

    public void setLoss(byte loss) {
        this.loss = loss;
    }

    public int getIonCode() {
        int lossOther = 5;
        return (((((int) type) * 20 + (int) charge)) * (lossOther + 1) + (int) loss) * 100 + (int) index + 1;
    }

    public static Product read(byte[] buffer, InputStream in) throws IOException {
        int bytesRead = in.read(buffer, 0, 12);
        if (bytesRead != 12) {
            throw new IOException("Unexpected end of stream while reading Product");
        }

        Product product = new Product();
        
        int mzBits = ((buffer[0] & 0xFF)) |
                     ((buffer[1] & 0xFF) << 8) |
                     ((buffer[2] & 0xFF) << 16) |
                     ((buffer[3] & 0xFF) << 24);
        product.mz = Float.intBitsToFloat(mzBits);
        
        int heightBits = ((buffer[4] & 0xFF)) |
                         ((buffer[5] & 0xFF) << 8) |
                         ((buffer[6] & 0xFF) << 16) |
                         ((buffer[7] & 0xFF) << 24);
        product.height = Float.intBitsToFloat(heightBits);
        
        product.charge = buffer[8];
        product.type = buffer[9];
        product.index = buffer[10];
        product.loss = buffer[11];
        
        return product;
    }
}

