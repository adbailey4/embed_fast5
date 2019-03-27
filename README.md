# embed_fast5
Embed fast5 signal files with a correctly formatted kmer to event alignment given the nucleotide sequence and a model file

####Install instructions  
* `git clone https://github.com/adbailey4/embed_fast5.git` 
* `cd embed_fast5`
* `python setup.py install`
* `mkdir cmake-build`
* `cd cmake-build`
* `cmake -DNANOPOLISH_HOME:STRING=/path/to/nanopolish ../`
* `make`
* `make check`

#### Usage
* `python run_embed_fast5.py --fast5 path/to/multi_read.fast5 --output_dir /path/to/output_dir --main_cpp_dir path/to/cmake-build/main_cpp --nanopolish_dir path/to/nanopolish --jobs 4`
