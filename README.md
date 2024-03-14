# DEploid_Pv_pipeline
Software Intellectual property belongs to: https://github.com/DEploid-dev/DEploid

This repo is just a user-friendly wrapper for analysing P. vivax data using DEploid software.

## Installation on local devices and private servers (skip this step if using ADA)
First create a tools directory to keep track of all software and repositories used by this pipeline.

    mkdir tools
    git clone https://github.com/pathogenseq/fastq2matrix.git
    cd fastq2matrix
    python setup.py install

    cd ..
    git clone https://github.com/aosborne13/DEploid_Pv_pipeline.git
    cd DEploid_Pv_pipeline
    python setup.py install

