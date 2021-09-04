[![Build Status](https://travis-ci.org/adbailey4/embed_fast5.svg?branch=master)](https://github.com/adbailey4/embed_fast5)

# embed_fast5

Embed fast5 signal files with a correctly formatted kmer to event alignment given the nucleotide sequence and a model file

#### Dependencies

* Boost >= 1_69_0 
* HDF5 (`sudo apt-get install libhdf5-dev`)
* Eigen (`sudo apt-get install libeigen3-dev`)
* HTSLIB >= 1.9
* cmake >= 3.15

### Install instructions  

* `git clone --recursive https://github.com/adbailey4/embed_fast5.git` 
* `cd embed_fast5`
* `pip install -e .`
* `pytest`


### Debugging

If the install instructions do not work, check out this [Dockerfile](https://github.com/adbailey4/base_images/blob/master/embed_dependencies/Dockerfile)
