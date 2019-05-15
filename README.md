[![Build Status](https://travis-ci.org/adbailey4/embed_fast5.svg?branch=master)](https://github.com/adbailey4/embed_fast5)


# embed_fast5
Embed fast5 signal files with a correctly formatted kmer to event alignment given the nucleotide sequence and a model file

#### Dependencies
* Boost (`sudo apt-get install libboost-all-dev`)
* HDF5 (`sudo apt-get install libhdf5-serial-dev`)
* Eigen (`sudo apt-get install libeigen3-dev`)
* OpenSSL (`sudo apt-get install libssl-dev`)
### Install instructions  


* `git clone https://github.com/adbailey4/embed_fast5.git` 
* `cd embed_fast5`
* `python setup.py install`
* `mkdir cmake-build`
* `cd cmake-build`
* `cmake ../`
* `make`
* `make check`

##### split_multi_fast5
split_multi_fast5 splits a multi_fast5_file into multiple fast5s. If there is a fastq object in the read, it will include that in the new individual files. split_multi_fast5 should be in your PATH after the install. 

* `split_multi_fast5 --multi_fast5_dir path/to/multi_read.fast5 --output_dir /path/to/output_dir --jobs 4`

|  Arguments | Help  | 
|---|---|
| -h, --help  | show this help message and exit  |
| --multi_fast5_dir MULTI_FAST5_DIR, -f MULTI_FAST5_DIR  |  path to multi read fast5 file |
| --output_dir OUTPUT_DIR, -o OUTPUT_DIR | Directory to place all new originally formatted fast5 files  |
| --jobs NB_JOBS, -j NB_JOBS |  number of jobs to run in parallel |
| --debug, -d  |  Option to not multiprocess reads to check for errors |
| --delete_multi | Option to delete the multi_fast5 file after processing  |


##### embed_fast5s
embed_fast5s will write a fastq and event table into a read using nanopolish's load from raw method. 
There is a preprocessing step which makes sure the fastq file is appropriately 
formatted and if it is not, will edit and write the new fastq file to the output directory.

* `embed_fast5s --fast5_dir path/to/fast5_dir  --fastq path/to/fastq/ --output_dir /path/to/output_dir --jobs 4 --embed_build_dir path/to/build_dir`

|  Arguments | Help  | 
|---|---|
| -h, --help  | show this help message and exit  |
| --fast5_dir FAST5_DIR, -f FAST5_DIR  |  path to signal fast5 files |
| --fastq FASTQ, -q FASTQ | path to fastq file  |
| --jobs NB_JOBS, -j NB_JOBS |  number of jobs to run in parallel |
| --debug, -d  |  Option to not multiprocess reads to check for errors |
| --output_dir OUTPUT_DIR, -o OUTPUT_DIR | Output directory for edited fastq if it is in wrong format|
| --embed_build_dir EMBED_BUILD_DIR, -m EMBED_BUILD_DIR  |  Directory where the embed_main executable is located |
| --no_fastq  |  If set this will not embed fastqs into the files |
| --no_events  |  If set this will not embed event tables into the files |


