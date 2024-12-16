# GeneGorman: Short Read Mapper

## Package Description
Read mapping is the process of aligning short DNA sequences (reads) to a reference genome to determine where they originate. This process is crucial for understanding genetic variations, which can reveal insights into human health, disease, and evolutionary biology. Short read mappers need to efficiently find matches for small DNA sequences within a much larger reference, even when there are many similar or repetitive regions.
To solve these problems, we created GeneGorman, an ode to the dean of students at the university of its creation, Rice. GeneGorman is a tool designed to address the challenges of short read mapping by providing a cost-effective, time-efficient, and memory-efficient solution. Named in honor of the dean of students at its founding university, Rice, GeneGorman parses input data, indexes the reference genome, and aligns reads using a seeding strategy. The system efficiently searches for exact matches of read segments (seeds) within the genome, making it possible to handle large-scale datasets with precision and speed.

## Team Members
- Delaney Miller
- Gal Kadmon
- Henry Tran
- Lauren Hu
- Yanni (Michelle) Pang

## Table of Contents
1. [Repository Structure](#repository-structure)
2. [Package Installation](#package-installation)
3. [Description](#description)
4. [Usage](#usage)
5. [Data](#data)
6. [Software Design](#software-design)
7. [Limitations/Known Issues](#limitations-and-known-issues)
8. [References](#references)

## Repository Structure
```
GeneGorman/
├── data/
│   ├── ground_truth.txt
│   ├── short_reads_1.fastq
│   ├── short_reads_2.fastq
│   └── short_reads_ref.fasta
├── preprocessing/
│   └── preprocesser.py
├── read_mapper/
│   └── read_mapper.py
├── output/
│   └── sam_output.py
├── validation/
│   └── validator.py
├── test/
│   ├── preprocessing_tests.py
│   └── read_mapper_tests.py
├── main.py
└── README.md
```

* `GeneGorman/` is the main directory for the GeneGorman genomic read mapper code.
    * `data/` contains the various data files used to demo and test GeneGorman.
	* `preprocessing/` includes code for parsing and preprocessing genomic data.
		* `preprocesser.py` – Code for preprocessing reads before mapping.
	* `read_mapper/` contains the mapper module that maps the reads to the reference genome.
		* `read_mapper.py` – Main Python code for the read mapping process
    * `output/` contains code for producing the sam output file after read mapping
        * `sam_output.py` - code for creating and outputting the sam file
	* `validation/` contains code for validating mapping results.
		* `validation.py` – code for validating and evaluating the accuracy of our resultant mapped reads.
    * `test/` contains code for testing individual modules
        * `preprocessing_tests.py` - code for testing the preprocessor
        * `read_mapper_tests.py` - code for testing the read mapper
	
## Package Installation

1. **Clone the Repository**: Within your VSCode environment clone the GitHub repository.

`$ git clone https://github.com/COMP-416-519-2024/GeneGorman.git`

2. **Navigate to the Project Directory**

`$ cd GeneGorman`

3. **Install Installation Prerequisites**: Use `pip` to install the BioPython

`pip install bio`

If you need to install under a specific version of Python, try something like this:

`python3.9 -m pip install biopython`
`pypy -m pip install biopython`

4. **Verify Installation**

    Open up a Python shell and try importing the package. If there are no errors, the package was successfully installed!    

```python

import GeneGorman

```

## Data

The instructions to acquire the dataset used to demo our tool can be found in `/data`.

### short_reads_ref.fasta

This file contains the reference genome sequence in FASTA format. It serves as the target against which the short reads are aligned.

### short_reads_1.fastq

This file contains the first set of short read sequences from paired-end sequencing. The FASTQ format includes both the nucleotide sequence and corresponding quality scores for each nucleotide.

### short_reads_2.fastq

Complementary to short_reads_1.fastq, this file holds the second set of short read sequences from the paired-end sequencing process. It also uses the FASTQ format, including both the sequence and quality score.

### ground_truth.txt

This file provides the known, correct alignment positions of certain reads, serving as a benchmark for validation. It is used to evaluate the accuracy of the mapping process by comparing the predicted alignments with these ground truth results.

## Software Design

### Program Input Data

The initial input consists of the raw genomic data: short reads and a reference genome. These data files must be in their standard formats—FASTQ for reads and FASTA for the reference genome. GeneGorman uses this data for preprocessing and alignment.

### Preprocessing and Indexing

The system parses and preprocesses the reads, extracting relevant data for efficient mapping. It then constructs a suffix array from the reference genome, enabling quick lookup of exact matches during the mapping phase.

### Seeding and Mapping

GeneGorman operates a "seed-and-extend" strategy to align reads. The mapper splits each read into smaller segments (seeds), searches for exact matches of these seeds in the genome using the suffix array, and identifies potential alignment positions. This strategy allows for faster, more accurate mapping by narrowing down the search space.

### Validation

GeneGorman assesses the performance of the mapping process by comparing results against a ground truth dataset. This helps ensure the accuracy and reliability of the system, providing metrics such as precision and recall to evaluate overall performance.

### Output

The final alignment results are compiled into a SAM (Sequence Alignment/Map) file. This standardized format ensures compatibility with other bioinformatics tools for further analysis.

## Limitations and Known Issues

GeneGorman’s current limitations stem from its use of the suffix array. Although efficient for most datasets, constructing the suffix array can be computationally intensive, especially for larger genomes. Additionally, the seed-based approach may struggle with highly repetitive regions, which can lead to ambiguous mappings. Future updates aim to address these limitations by exploring optimization strategies and alternative indexing methods.

## References

Bowe, Alex. “FM-Indexes and Backwards Search.” Alexbowe.Com, alexbowe.com, 22 Feb. 2023, www.alexbowe.com/fm-index/. 

Musich R, Cadle-Davidson L, Osier MV. Comparison of Short-Read Sequence Aligners Indicates Strengths and Weaknesses for Biologists to Consider. Front Plant Sci. 2021;12: 657240.

Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9: 357–359.

Nate. “Short Read Alignment: Seeding - Seven Bridges Genomics.” Seven Bridges, Seven 
Bridges, 4 Oct. 2018, www.sevenbridges.com/short-read-alignment-seeding/#:~:text=Alignment%3A%20a%20quick%20review&text=Many%20modern%20alignment%20algorithms%20rely,so%20seeding%20is%20very%20fast. 

Langmead, Ben. “Langmead Lab @ JHU - Teaching Materials.” Langmead-Lab.org, 2023, www.langmead-lab.org/teaching.html.

N. Ahmed, K. Bertels and Z. Al-Ars, "A comparison of seed-and-extend techniques in modern DNA read alignment algorithms," 2016 IEEE International Conference on Bioinformatics and Biomedicine (BIBM), Shenzhen, China, 2016, pp. 1421-1428, doi: 10.1109/BIBM.2016.7822731. keywords: {DNA;Computers;Indexes;Genomics;Bioinformatics;Irrigation},