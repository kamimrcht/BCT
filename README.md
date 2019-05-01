# BCT
## de Bruijn graph based corrector for transcriptomic data

## Installation:
```./install.sh -t 8```

Will install bct using 8 thread for the compilation

## Usage:

Basic usage:

```./Bct.py -u reads.fa -o working_directory -t core_number ```

More agressive filtering:

```./Bct.py -u reads.fa -o working_directory -t core_number  -S abundance_threshold```

With this option the graph will contain less errors but regions seen less than S times may be lost

Disable poly A tail cleaning:

```./Bct.py -u reads.fa -o working_directory -t core_number  -S abundance_threshold -c 0```

Use this option if you do not have polyA tails in your dataset, in meta-genomic for example.

The polyA tail is HIGHLY recommanded for transriptomic data.

Other options are advanced parameters, change them at your own risk


## Option:

#### -h display help 
You are lost

#### -u read file 
Your input reads in fasta/multifasta/fastq file gziped or not

If you have multiple dataset please concatenate them

I assume fastq file will contain ".fq" or ".FQ" or ".fastq" or ".FASTQ" if not please add ```-q``` to the command line

#### -o output directory
The cleaned graph, the logs and the corrected reads will be put in this directory (default .)

#### -t thread number
The number of threads we can run in parallel (default 20)

#### -k kmer size
The kmer size used to build the graph (default 31)

#### -s kmer abundance threshold
kmer seen less than s time will be lost (default 2)

#### -S unitig abundance threshold
unitig seen less than S time will be lost (default 2)

#### -i anchors subsampling
 Index 1 out of i anchors to reduce memory usage (default all)
 
 #### -n repeated anchors
 Anchors seen more than n time are not indexed (default 8) 
 
 #### -d debug mode
 Will print the command line used by BCT
 
 
 


