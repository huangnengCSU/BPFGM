## BlockPolish feature generation module
BlockPolish is an accurate polishing of long-read assembly via block divide-and-conquer. In the method, BlockPolish divides contigs into blocks with low complexity and high complexity according to statistics of reads aligned to the assembly. Multiple sequence alignment is applied to re-align raw reads in complex blocks and improve the consistency of the alignment result. Due to the different distributions of error rates in trivial and complex blocks, two multi-task bi-directional LSTM networks are proposed to predict the consensus sequences.

BlockPolish feature generation module is used to divide the contigs into trivial blocks and complex blocks from reads-to-assembly alignment ```BAM``` file. Then, this module generates a feature matrix for each block by counting the base distribution at each position. The features of trivial blocks and complex blocks are written to two output file named ```TRIVIAL_BLOCK_FEATURES``` and  ```COMPLEX_BLOCK_FEATURES```.

## Installation
---
### Dependencies
if compiling on Ubuntu, this will install all required packages:
```
apt-get -y install make gcc g++ zlib1g-dev
```

BlockPolish feature generation module is compiled with cmake. We recommend using the latest cmake version, but 3.10 and higher are supported:
```
wget https://github.com/Kitware/CMake/releases/download/v3.19.4/cmake-3.19.4-Linux-x86_64.sh && sudo mkdir /opt/cmake && sudo sh cmake-3.19.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && sudo ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
cmake --version
```

### Compilation
- Check out repository:
```
git clone https://github.com/huangnengCSU/BPFGM.git 
cd BPFGM
```

- Make build directory:
```
mkdir build
cd build
```

- Generate Makefile and run:
```
cmake ..
make
./block
```

## Running BlockPolish feature generation module
```
usage: block -b <BAM_FILE> -s <TRIVIAL_BLOCK_FEATURES> -c <COMPLEX_BLOCK_FEATURES>
Generate the features for BlockPolish polishing using reads-to-assembly alignments in BAM_FILE.

Required arguments:
    BAM_FILE is the alignment of reads to the assembly.
    TRIVIAL_BLOCK_FEATURES is the output file of the features of the trivial blocks.
    COMPLEX_BLOCK_FEATURES is the output file of the features of the trivial blocks.
```

### Sample Execution
```
./block -b ../test/HG002.chr22.2k.bam -s HG002.chr22.2k.trivial_features.txt -c HG002.chr22.2k.complex_features.txt
```

## Data Formats    

BlockPolish feature generation module requires that the input BAM is indexed.

---
## License
Copyright (C) 2020 by Huangneng (huangn94@foxmail.com)
