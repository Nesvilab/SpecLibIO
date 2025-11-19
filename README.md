# SpecLibIO

A Java library for reading and writing DIA-NN spectral library files (`.speclib` format).

## Overview

SpecLibIO provides a pure Java implementation for working with DIA-NN spectral libraries, enabling seamless integration with Java-based proteomics tools and workflows. The library supports reading and writing the binary `.speclib` format used by [DIA-NN](https://github.com/vdemichev/DiaNN).

## Development

This project, **including all implementation code, unit tests, and this README**, was generated using **vibe coding** (AI-assisted development). The code has been reviewed, (partially) modified, and tested by Fengchao Yu.

## Attribution

### Reader Implementation
The `DiaNNSpecLibReader` is based on the DIA-NN SpecLib parser from the [ProteoWizard](https://github.com/ProteoWizard/pwiz) project.

### Writer Implementation
The `DiaNNSpecLibWriter` is based on the implemented `DiaNNSpecLibReader`, providing symmetric write capabilities for the binary format.

### Test Data
Unit test files (`.speclib` and `.tsv` files) are from the [ProteoWizard](https://github.com/ProteoWizard/pwiz) repository's test resources.

## Features

- **Read DIA-NN spectral libraries**: Parse binary `.speclib` files into Java objects
- **Write DIA-NN spectral libraries**: Generate binary `.speclib` files from Java objects
- **Parquet to SpecLib conversion**: Convert Parquet-formatted spectral libraries to SpecLib format
- **Version support**: Supports format versions up to -3 (latest DIA-NN format)
- **Complete data model**: Includes proteins, peptides, precursors, fragments, and metadata
- **Stream-based I/O**: Efficient handling of large library files
- **Pure Java**: No native dependencies, works on any platform with Java 11+

## Requirements

- Java 11 or higher
- Maven 3.x (for building)

## Usage

### Reading a Spectral Library

```java
import speclib.io.DiaNNSpecLibReader;
import speclib.io.SpectralLibrary;

// From file path
DiaNNSpecLibReader reader = new DiaNNSpecLibReader("library.speclib");
SpectralLibrary library = reader.read();

// From input stream
try (InputStream is = new FileInputStream("library.speclib")) {
    DiaNNSpecLibReader reader = new DiaNNSpecLibReader(is);
    SpectralLibrary library = reader.read();
}
```

### Writing a Spectral Library

```java
import speclib.io.DiaNNSpecLibWriter;
import speclib.io.SpectralLibrary;

SpectralLibrary library = new SpectralLibrary();
library.setName("My Library");
library.setGenDecoys(true);
library.setGenCharges(true);
// ... populate library with entries

DiaNNSpecLibWriter writer = new DiaNNSpecLibWriter(library);
writer.write("output.speclib");
```

### Converting from Parquet

```java
import speclib.io.ParquetToSpecLib;

// Convert Parquet file to SpecLib
ParquetToSpecLib converter = new ParquetToSpecLib("input.parquet");
converter.convertAndWrite("output.speclib");
```

## Testing

Run unit tests with Maven:

```bash
mvn test
```

The test suite includes:
- **DiaNNSpecLibReaderTest**: Reading various DIA-NN library formats
- **DiaNNSpecLibWriterTest**: Round-trip read/write validation, edge cases, format version compatibility
- **ParquetToSpecLibTest**: Parquet to SpecLib conversion with comprehensive field validation
  - Tests conversion accuracy for all fields (precursors, fragments, proteins, RT, ion mobility)
  - Validates fragment types, loss types, charges, and intensities
  - Verifies round-trip conversion (Parquet → SpecLib → Read back)

## Binary Format Details

The DIA-NN `.speclib` format is a binary format with the following structure:

1. **Header**: Version, flags (gen_decoys, gen_charges, infer_proteotypicity)
2. **Metadata**: Library name, FASTA names
3. **Proteins**: Array of protein isoforms and protein groups
4. **Precursors**: Array of precursor identifiers
5. **Names/Genes**: Parallel arrays for protein names and gene names
6. **iRT range**: Minimum and maximum indexed retention times
7. **Entries**: Array of library entries with peptides and fragments
8. **Optional**: Elution groups (version -1 and later)

All integers are stored as 4-byte little-endian, doubles as 8-byte little-endian IEEE 754.

