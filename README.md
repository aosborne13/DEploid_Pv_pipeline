# DEploid_Pv_pipeline
Software Intellectual property belongs to: https://github.com/DEploid-dev/DEploid

This repo is just a user-friendly wrapper for analysing P. vivax data using DEploid software.

## Installation on local devices and private servers
First create a tools directory to keep track of all software and repositories used by this pipeline.

    mkdir tools
    cd tools

Clone required GitHub Repositories

    git clone https://github.com/pathogenseq/fastq2matrix.git
    cd fastq2matrix
    python setup.py install

    cd ..
    git clone https://github.com/aosborne13/DEploid_Pv_pipeline.git
    cd DEploid_Pv_pipeline
    python setup.py install

