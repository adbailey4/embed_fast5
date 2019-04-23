[![Build Status](https://travis-ci.org/adbailey4/embed_fast5.svg?branch=master)](https://github.com/adbailey4/embed_fast5)


# embed_fast5
Embed fast5 signal files with a correctly formatted kmer to event alignment given the nucleotide sequence and a model file

####Install instructions  
* `git clone https://github.com/adbailey4/embed_fast5.git` 
* `cd embed_fast5`
* `python setup.py install`
* `mkdir cmake-build`
* `cd cmake-build`
* `cmake ../`
* `make`
* `make check`

#### Usage
* `python run_embed_fast5.py --fast5 path/to/multi_read.fast5 --output_dir /path/to/output_dir --main_cpp_dir path/to/cmake-build/main_cpp --nanopolish_dir path/to/nanopolish --jobs 4`

##### split_multi_fast5
split_multi_fast5 splits a multi_fast5_file into multiple fast5s. If there is a fastq object in the read, it will include that in the new individual files. split_multi_fast5 should be in your PATH after the install. 

* `split_multi_fast5 --fast5_dir path/to/multi_read.fast5 --output_dir /path/to/output_dir --jobs 4`

|  Arguments | Help  | 
|---|---|
| -h, --help  | show this help message and exit  |
| --fast5_dir FAST5_DIR, -f FAST5_DIR  |  path to multi read fast5 file |
| --output_dir OUTPUT_DIR, -o OUTPUT_DIR | Directory to place all new originally formatted fast5 files  |
| --jobs NB_JOBS, -j NB_JOBS |  number of jobs to run in parallel |
| --debug, -d  |  Option to not multiprocess reads to check for errors |
| --delete_multi | Option to delete the multi_fast5 file after processing  |
