#  rrequested-pkg 

### -General purpose and applicability

rrequested-pkg is the python package distribution of the modular shellscript tool RREQUESTED (Raw REads QUality Extraction, Size Trimming and Ex-novo Demultiplexing). It is aimed at preprocessing (quality filtering, size selection and demultiplexing) of raw basecalled reads, especially produced by Oxford Nanopore, a Third Generation Sequencing technology.

Even though the program has been tested only on Nanopore sequencing results, the program could be cautiously applied also to products from NGS or other ThirdGen technologies. 


### -Installation 
The basic requirment for this installation is to have python 3.10 installed on your machine.

To install the package, just run:
```
pip install rrequested-pkg
```

You can retrieve general information on it by running:
```
pip show rrequested-pkg
```

### -Functions and usage
| Module | Class | Method | Description |
| --------- | --------- | --------- | --------- |
| **file_system** | _FileSystem_ | `makedir_orchange` | Method returns the provided path after attempting to create the directory, handling the case where the directory already exists. |
| **file_system** | _FileSystem_ | `get_base_dir` | Method extracts information such as the base directory, base name, and the base directory of the base directory from the provided file path, using the os.path.splitext function and string manipulation. |
| **quality_filter** | _Quality_Filter_ | `load_data` | Method reads data from the input file, handling both gzipped and non-gzipped fastq formats. |
| **quality_filter** | _Quality_Filter_ | `ascii_conv_and_mean` | Method converts ASCII characters to Phred quality scores and calculates the mean. |
| **quality_filter** | _Quality_Filter_ | `filter` | Method performs the filtering based on the quality threshold and updates the input file with filtered reads, providing a summary of the process. |
| **size_filter** | _Size_Filter_ | `load_data` | Method reads data from the input file, handling both gzipped and non-gzipped fasta/fastq formats. |
| **size_filter** | _Size_Filter_ | `size_filter` | Method performs the size-based filtering and writes the selected sequences to a new file, providing a summary of the process. |
| **unref_demultiplex** | _Unreferenced_demultiplexing_ | `find_the_num` | Method finds highly divergent sequences from the input file. |
| **unref_demultiplex** | _Unreferenced_demultiplexing_ | `find_the_num_list` | Method finds highly divergent sequences from a given list of sequences. |
| **unref_demultiplex** | _Unreferenced_demultiplexing_ | `demultiplex` | Method performs the demultiplexing process and writes the demultiplexed sequences into separate files. |

Here are some example usages of the package:

**Size_Filter**

```python
# Import necessary modules
from rrequested_pkg.size_filter import Size_Filter

# Specify input file, minimum size, and maximum size
input_file = "your_input_file.fastq"
min_size = 20
max_size = 200

# Create an instance of Size_Filter
size_filter_instance = Size_Filter(infile=input_file, minsize=min_size, maxsize=max_size)

# Load data from the input file
data = size_filter_instance.load_data()

# Perform size filtering
filtered_data = size_filter_instance.size_filter()
```

**Quality_Filter**

```python
# Import necessary modules
from rrequested_pkg.quality_filter import Quality_Filter

# Specify input file and quality threshold
input_file = "your_input_file.fastq"
quality_threshold = 20

# Create an instance of Quality_Filter
quality_filter_instance = Quality_Filter(infile=input_file, quality=quality_threshold)

# Load data from the input file
data = quality_filter_instance.load_data()

# Perform quality filtering
quality_filter_instance.filter()
```

**FileSystem**
```python
# Import necessary modules
from rrequested_pkg.file_system import FileSystem

# Create an instance of FileSystem
file_system_instance = FileSystem()

# Specify the path where you want to create a directory
directory_path = "/path/to/your/directory"

# Create a new directory or change to an existing one
base_directory = file_system_instance.makedir_orchange(directory_path)

# Get information about the base directory
base_dir, base_name, base_basedir = file_system_instance.get_base_dir(base_directory)

# Display the directory information
print("Base Directory:", base_dir)
print("Base Name (without extension):", base_name)
print("Base Directory of the Base Directory:", base_basedir)
```

**Unreferenced_Demultiplexing**
```python
# Import necessary modules
from rrequested_pkg.unreferenced_demultiplexing import Unreferenced_Demultiplexing

# Specify input file
input_file = "input_file.fasta"

# Create an instance of Unreferenced_Demultiplexing
demultiplex_instance = Unreferenced_Demultiplexing(infile=input_file)

# Find highly divergent sequences
highly_divergent_seqs = demultiplex_instance.find_the_num()

# Demultiplex the sequences
demultiplex_instance.demultiplex()

# The demultiplexed sequences will be written to individual files in a new folder
# The folder will be named based on the input file, in this example: "input_file-fa-demultiplexed"
```

### -How does it work? ###
1. The quality filtering method is based on the easiest implementation one could think of: for every read, the filtering algorithms takes the mean quality and discards the reads that are under a given value (default is 7, so this step will take place nevertheless if the file is fastq/fastq.gz)
4. The size filtering method is also based on the easiest implementation one could think of: for every read, the filtering algorithms takes the length and, if this is below the minimum or above the maximum allowed, the read gets discarded. 
5. The demultiplexing method is based on super-fast global alignment and it is divided into two main parts: in the first, the demultiplexer identifies unique (higly divergent) reads, that are the ones which score less than 50% in similarity with all the other sequences. After that, it globaly aligns all the raw reads against the "self-made reference", grouping the ones that share more than 70% of their code. After having demultiplexed this way, the program checks the leftovers, to see wether there are worthy-to-save data or not. This brings to a five-round cycle that identifies higly divergent sequences in the non-grouped ones and clusters the "nogroup" reads against them. If there are still ungrouped remainders after this step, they get clustered together (3). Only groups encompassing more than 1% of the total reads will be written as demultiplexed fasta files, named N.fasta (where N is a number) or nogroup.fasta if they belong to the unclustered reads: they could be found in the folder basefilename-extensionabbreviation-demultiplexed (an example could be: if you are demultiplexing a file named coleoptera.fastq, the folder will be coleoptera-fq-demultiplexed)



### -Final considerations ###
As a practical suggestion, we strongly advise to be cautious while using rrequested-pkg with files containing reads from multiple individuals, especially if the quality of the data is low: not because there is the risk that you will miss something, but because it can produce more groups than needed.

Moreover, please note that rrequested-pkg is still experimental and may contain errors, may fail/take really long while performing large analyses with limited computational power (e.g. on a normal laptop) and may output not-100%-reliable results, so always check them and pull issues whenever you feel it to be the case, we'll be on your back as soon as possible to fix/implement/enhance whatever you suggest!


### -License and rights of usage
The code is distributed under the MIT license.

The MIT License is a straightforward and permissive license that encourages collaboration and widespread use. It allows developers the freedom to build upon and share code while providing clear guidelines for responsible and respectful use. Whether you are an individual or a corporation, the MIT License promotes an open and collaborative approach to software development.