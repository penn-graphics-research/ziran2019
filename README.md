## Code description 

This is the opensource code for the following papers:

(1) Silly Rubber: An Implicit Material Point Method for Simulating Non-equilibrated Viscoelastic and Elastoplastic Solids ,Yu Fang, Minchen Li, Ming Gao, Chenfanfu Jiang, (SIGGRAPH 2019)

(2) CD-MPM: Continuum Damage Material Point Methods for Dynamic Fracture Animation ,Joshuah Wolper, Yu Fang, Minchen Li, Jiecong Lu, Ming Gao, Chenfanfu Jiang, (SIGGRAPH 2019)

It is tested on a fresh install of Ubuntu 18.04 LTS.

## Unzip Data

Go to Data/LevelSets and unzip **breadxxx.vdb.zip** into the same directory.

You need to do this due to the github single file size limit.

## Dependencies Installation

    sudo apt-get install make cmake g++ libeigen3-dev gfortran libmetis-dev
    sudo apt-get install libopenvdb-dev libboost-all-dev libilmbase-dev libopenexr-dev
    sudo apt-get install libtbb2 libtbb-dev libz-dev clang-format-6.0 clang-format
   
## Building in Ziran

    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j 4

## Running Demos

    cd Projects/mpm
    ./mpm -test 1
    ./mpm -test 2
    ./mpm -test 3

    cd Projects/admm
    ./admm -test 1
    ./admm -test 2
    ./admm -test 3
    ./admm -test 4
    ./admm -test 5
    ./admm -test 6
    
    cd Projects/fracture
    ./fracture -test 1
    ./fracture -test 2
    ./fracture -test 3

## Bibtex

Please cite our papers if you use this code for your research: 
```
@article{fang2019silly,
  title={Silly rubber: an implicit material point method for simulating non-equilibrated viscoelastic and elastoplastic solids},
  author={Fang, Yu and Li, Minchen and Gao, Ming and Jiang, Chenfanfu},
  journal={ACM Transactions on Graphics (TOG)},
  volume={38},
  number={4},
  pages={118},
  year={2019},
  publisher={ACM}
}
```
```
@article{wolper2019cd,
  title={CD-MPM: Continuum damage material point methods for dynamic fracture animation},
  author={Wolper, Joshuah and Fang, Yu and Li, Minchen and Lu, Jiecong and Gao, Ming and Jiang, Chenfanfu},
  journal={ACM Transactions on Graphics (TOG)},
  volume={38},
  number={4},
  pages={119},
  year={2019},
  publisher={ACM}
}
```
