<p align="center">
   <img src="https://github.com/vinni-bio/EXONtools/blob/master/img/EXONtools_small.png" alt="EXONtoools logo" width="600" height="115">
</p>

<h2 align="center">version 0.3b</h2>

<p align="center">
  <strong>A complete and flexible pipeline for exon capture sequencing data analysis of non-model organisms</strong>
  <br>
  <a href="https://raw.githubusercontent.com/vinni-bio/EXONtools/master/tutorial/EXONtools_v01b.pdf">Download Tutorial</a>
  ·
  <a href="https://github.com/vinni-bio/EXONtools/blob/master/readme/pipeline-scheme.md">Pipeline</a>
  ·
  <a href="https://github.com/vinni-bio/EXONtools/blob/master/readme/list-of-commands.md">Commands</a>
  ·
  <a href="https://github.com/vinni-bio/EXONtools/issues/new">Report bug</a>
</p>


## Table of contents

- [Exon capture sequencing](https://github.com/vinni-bio/EXONtools#exon-capture-sequencing)
- [EXONtools project](https://github.com/vinni-bio/EXONtools#exontools-project)
- [Quick start](https://github.com/vinni-bio/EXONtools#quick-start)
- [INSTALLATION](https://github.com/vinni-bio/EXONtools#installation)
- [Pipeline scheme](https://github.com/vinni-bio/EXONtools#pipeline-scheme)
- [EXONtools commands](https://github.com/vinni-bio/EXONtools#exontools-commands-dependencies)
- [EXONtools general options](https://github.com/vinni-bio/EXONtools#exontools-general-options)
- [Supplementary EXONtools commands](https://github.com/vinni-bio/EXONtools#supplementary-exontools-commands)
- [Links to EXONtools dependencies](https://github.com/vinni-bio/EXONtools#links-to-exontools-dependencies)
- [Step-by-step example](https://github.com/vinni-bio/EXONtools#step-by-step-example)
- [Future development](https://github.com/vinni-bio/EXONtools#future-development)
- [Glossary](https://github.com/vinni-bio/EXONtools#glossary)
- [Citation](https://github.com/vinni-bio/EXONtools#Citation)
- [Acknowledgements](https://github.com/vinni-bio/EXONtools#acknowledgements)
- [Contacts](https://github.com/vinni-bio/EXONtools#econtacts)


## Exon capture sequencing

Exon capture is one of the most promising approaches in next generation sequencing. By targeting only protein coding regions it saves a lot of sequencing effort and allows researchers to obtain SNP data for thousands of genes in hundreds of individuals during one high-throughput sequencing run. Exon capture sequencing is especially important for evolutionary and population genomic studies of non-model organisms (Genotyping by Sequencing) because:
- this method does not require a reference genome
- this method does not rely solely on conservative regions like UCE
- this method allows to obtain SNPs from both exon (often non-neutral) and intron (often neutral) regions
- this method provides sufficient amount of SNP data to deep phylogenomic reconstructions and large scale population studies
- this method allows to screen tens of thousands of loci (genes) to identify those that are under selection

<p align="center">
   <img src="https://github.com/vinni-bio/EXONtools/blob/master/img/ECS_small.png" alt="Exon Capture Sequencing" width="800" height="250">
   <i><b>Figure 1.</b> Three major stages of laboratory protocol for exon capture sequencing project on a non-model organism.</i>
</p>

## EXONtools project

The main goal of the EXONtools pipeline is to help researchers who are working on their own exon capture sequencing (ECS) project to:
- reconstruct the complete annotated transcriptome of any non-model organism to create a "pseudoreference"
- prepare a set of biotin-labelled oligonucleotide sequences (hybridization baits) that will be used for target DNA enrichment
- produce the final dataset of SNPs based on a comprehensive analysis of ECS data collected from thousands of genes and hundreds of individuals

Written in Python, the EXONtools pipeline combines many well-known bioinformatics programs (dependencies) that are commonly used in standard NGS analyses and also implements many new algorithms within a single analytical framework.

## Quick start

Please check that the names of your sequence files follow the EXONtools rules of library name format:

**unpaired files:** LIBNAME_anyletters.fq OR LIBNAME.fasta (e.g., `ST-1_unpaired.fq`)<br/>
**paired files:** LIBNAME_R1.fq OR LIBNAME_R2.fq (e.g., `ST-1_R1_paired.fq`, `ST-1_R2_paired.fq`)

`R1` and `R2` annotations are mandatory for file names with paired reads. Otherwise those files be considered as unpaired. Each library (e.g., sample) can only have one file
with unpaired reads OR two files with paired reads (`R1` & `R2`).

Every computational task performed in the EXONtools pipeline is considered as an isolated step and therefore must be started with a separate command. Such an approach provides the end user with a flexibility in the adjustment of different pipeline parameters for each ECS data analysis, including changes in task order or omitting any unnecessary tasks. The EXONtools program is also designed to achieve the most efficient use of all available hardware resources. Each step of the EXONtools pipeline is supplemented by a comprehensive summary report that allows the user to compare and to optimize produced results.

**Each EXONtools command is executed using the following syntax:**
```
EXONtools.py [general options] COMMAND [command options]
```

The [general options](https://github.com/vinni-bio/EXONtools#exontools-general-options) of the EXONtools pipeline regulate system memory usage, multithreading, debugging modes and the output verbosity, including summary reports with comprehensive statistical information about each procedure. Additionally, the user can adjust some dependency program specific parameters by appending its command line arguments in the EXONtools general options (`--extra`).

Each EXONtools command has a unique set of its own arguments ([command options](https://github.com/vinni-bio/EXONtools#exontools-commands-dependencies)). Please read them carefully before launching a task.

## Installation

The EXONtools pipeline does not require any special compilation in your system environment and can be directly cloned from the EXONtools repository:

```
git clone https://github.com/vinni-bio/EXONtools.git
```

or [download the latest release](https://github.com/vinni-bio/EXONtools/archive/v0.2b.zip)

Please install the following Python packages with PIP:
- [NUMPY](https://numpy.org/install) `pip install numpy`
- [PSUTIL](https://psutil.readthedocs.io/en/latest/) `pip install psutil`

**IMPORTANT!!!** Many steps in the EXONtools pipeline rely on one or several dependency programs that must be preinstalled manually ~~or by running the docker image of the EXONtools pipeline~~ (*in prep.*). The list of all required dependencies is provided [here](https://github.com/vinni-bio/EXONtools#links-to-exontools-dependencies) and in the configuration file `dependencies.ini` located in `EXONtools/src` directory. **Please provide the paths for all required dependency programs in the `dependencies.ini` file before using the pipeline.**

**IMPORTANT!!!** Installation of all dependencies is not required. Some dependencies can be omitted if the corresponding pipeline step is not going to be used in the analysis of ECS data. To skip dependency installation just leave a blank line in the \[PATHS\] option of the `dependencies.ini` file like that:

```fastqc =```

**IMPORTANT!!!** If dependency is installed in the environment PATH  (e.g., `/usr/local/bin`), just leave the `default` value in the \[PATHS\] option of the `dependencies.ini` file. In this case, the default command must **match the program name**:

```fastqc = default```

**IMPORTANT!!!** The EXONtools pipeline is currently designed to work only on Linux and Unix based operating systems with the command line interface. Supported Python language versions are 2.7.>10 or 3.5.>5.

## Pipeline scheme

```text
EXONtools.py
└── LAB PROTOCOL I: De novo transcriptome sequencing
        └── Stage A. Processing raw reads
               ├── Step A1. Demultiplexing (demultiplex_reads)
               ├── Step A2. Formatting (format_reads)
               ├── Step A3. Error correction (correct_reads)
               ├── Step A4. Duplication removal (deduplicate_reads)
               ├── Step A5. Trimming and filtering (clean_reads)
               ├── Step A6. Filtering low complexity (filter_reads)
               ├── Step A7. Merging paired reads (merge_reads)
               ├── Step A8. Contamination removal (decontaminate_reads)
               ├── Step A9. Formatting (format_reads)
               └── Step A10. Evaluating reads (evaluate_reads)

        └── Stage B. Pseudoreference construction
               ├── Step B1. Assembling reads (assemble_reads)
               ├── Step B2. Generating consensus assembly (consensus_assembly)
               ├── Step B3. Mapping reads to assembly (map_reads)
               ├── Step B4. Base calling (call_bases)
               ├── Step B5. Assembly annotation (annotate_contigs)
               ├── Step B6. Evaluate mapping quality (evaluate_mapping)
               └── Step B7. Evaluate assembly quality (evaluate_assembly)

         └── Stage C. Hybridization bait design
               ├── Step C1. Predicting exon boundaries (search_exons)
               ├── Step C2. Mapping exons to pseudoreference (map_exons)
               └── Step C3. Designing hybridization baits (design_baits)

└── LAB PROTOCOL II: Bait synthesis and exon capture

└── LAB PROTOCOL III: Exon capture sequencing
        └── Stage D. Processing raw reads
               ├── Step D1. Demultiplexing (demultiplex_reads)
               ├── Step D2. Formatting (format_reads)
               ├── Step D3. Error correction (correct_reads)
               ├── Step D4. Duplication removal (deduplicate_reads)
               ├── Step D5. Trimming and filtering (clean_reads)
               ├── Step D6. Filtering low complexity (filter_reads)
               ├── Step D7. Merging paired reads (merge_reads)
               ├── Step D8. Contamination removal (decontaminate_reads)
               ├── Step D9. Formatting (format_reads)
               └── Step D10. Evaluating reads (evaluate_reads)

        └── Stage E. SNP dataset construction
               ├── Step E1. Assembling reads (assemble_reads)
               ├── Step E2. Generating consensus assembly (consensus_assembly)
               ├── Step E3. Mapping reads to assembly (map_reads)
               ├── Step E4. Base calling (call_bases)
               ├── Step E5. Assembly annotation (annotate_contigs)
               ├── Step E6. Clustering exons (stack_exons)
               ├── Step E7. Aligning exons (align_stacks)
               ├── Step E8. Trimming and filtering exons (clean_stacks)
               ├── Step E9. Predicting exon boundaries (split_stacks)
               ├── Step E10. SNP calling from exon regions (call_snps)
               ├── Step E11. SNP calling from intron regions (call_snps)
               └── Step E12. SNP filtering (call_snps)
```

## [EXONtools commands](https://github.com/vinni-bio/EXONtools/blob/master/readme/list-of-commands.md): [dependencies](https://github.com/vinni-bio/EXONtools#links-to-exontools-dependencies)

1. [align_stacks](https://github.com/vinni-bio/EXONtools/blob/master/readme/align_stacks.md):
   - MAFFT ([Katoh and Standley, 2013](https://doi.org/10.1093/molbev/mst010))
2. [annotate_contigs](https://github.com/vinni-bio/EXONtools/blob/master/readme/annotate_contigs.md):
   - BLAST ([Altschul et al., 1990](https://doi.org/10.1016/S0022-2836(05)80360-2))
3. [assemble_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/assemble_reads.md): *choose one*
   - ABySS ([Simpson et al., 2009](https://doi.org/110.1101/gr.089532.108))
   - SPAdes ([Bankevich et al., 2012](https://doi.org/10.1089/cmb.2012.0021))
   - TransABySS ([Robertson et al., 2010](https://doi.org/10.1038/nmeth.1517))
   - Trinity ([Grabherr et al., 2011](https://doi.org/10.1038/nbt.1883))
4. [call_bases](https://github.com/vinni-bio/EXONtools/blob/master/readme/call_bases.md)
5. [call_snps](https://github.com/vinni-bio/EXONtools/blob/master/readme/call_snps.md)
6. [clean_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/clean_reads.md): *choose one*
   - BBduk (Bushnell, 2018)
   - Cutadapt ([Martin, 2011](https://doi.org/10.14806/ej.17.1.200))
   - Trimmomatic ([Bolger et al., 2014](https://doi.org/10.1093/bioinformatics/btu170))
7. [clean_stacks](https://github.com/vinni-bio/EXONtools/blob/master/readme/clean_stacks.md)
8. [consensus_assembly](https://github.com/vinni-bio/EXONtools/blob/master/readme/consensus_assembly.md):
   - BLAT ([Kent, 2002](https://doi.org/10.1101/gr.229202))
   - CAP3 ([Huang & Madan, 1999](https://doi.org/10.1101/gr.9.9.868))
   - CD-HIT ([Li & Godzik, 2006](https://doi.org/10.1093/bioinformatics/btl158))
9. [correct_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/correct_reads.md):
   - BayesHammer ([Nikolenko et al., 2013](https://doi.org/10.1186/1471-2164-14-S1-S7))
10. [decontaminate_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/decontaminate_reads.md):
    - Bowtie2 ([Langmead & Salzberg, 2012](https://doi.org/10.1038/nmeth.1923))
11. [deduplicate_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/deduplicate_reads.md)
12. [demultiplex_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/demultiplex_reads.md)
13. [design_baits](https://github.com/vinni-bio/EXONtools/blob/master/readme/design_baits.md)
14. [evaluate_assembly](https://github.com/vinni-bio/EXONtools/blob/master/readme/evaluate_assembly.md)
15. [evaluate_mapping](https://github.com/vinni-bio/EXONtools/blob/master/readme/evaluate_mapping.md)
16. [evaluate_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/evaluate_reads.md)
17. [filter_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/filter_reads.md)
18. [format_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/format_reads.md)
19. [map_exons](https://github.com/vinni-bio/EXONtools/blob/master/readme/map_exons.md):
    - BLAST ([Altschul et al., 1990](https://doi.org/10.1016/S0022-2836(05)80360-2))
    - CD-HIT ([Li & Godzik, 2006](https://doi.org/10.1093/bioinformatics/btl158))
20. [map_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/map_reads.md):
    - Bowtie2 ([Langmead & Salzberg, 2012](https://doi.org/10.1038/nmeth.1923))
    - bwa ([Li & Durbin, 2009](https://doi.org/10.1093/bioinformatics/btp324))
    - SAMtools ([Li et al., 2009](https://doi.org/10.1093/bioinformatics/btp352))
21. [merge_reads](https://github.com/vinni-bio/EXONtools/blob/master/readme/merge_reads.md): *choose one*
    - BBmerge ([Bushnell et al., 2017](https://doi.org/10.1371/journal.pone.0185056))
    - FLASH ([Magoč & Salzberg, 2011](https://doi.org/10.1093/bioinformatics/btr507))
22. [search_exons](https://github.com/vinni-bio/EXONtools/blob/master/readme/search_exons.md):
    - BLAST ([Altschul et al., 1990](https://doi.org/10.1016/S0022-2836(05)80360-2))
23. [split_stacks](https://github.com/vinni-bio/EXONtools/blob/master/readme/split_stacks.md)
24. [stack_exons](https://github.com/vinni-bio/EXONtools/blob/master/readme/stack_exons.md)

## EXONtools general options

- `-c`, `--citations`: show citations for EXONtools and for all dependencies
- `-D`, `--debug`: debugging mode for developers
- `-E`, `--extra`: auxillary command line arguments for dependency programs (leave a blank space before the opening quote, otherwise the EXONtools command will raise an error).<br/>Example: `EXONtools.py --extra " -argument value" COMMAND [command options]`. Recommended for advanced users only.
- `-K`, `--keeptmp`: do not delete temporary files after finishing task
- `-L`, `--log`: additionally save all console output in `EXONtools.log` file within the current directory
- `-l`, `--license`: show EXONtools license file
- `-M`, `--memory`: set maximum RAM usage for dependency programs (Gb)
- `-Q`, `--quiet`: hide console output except critical error messages (quiet mode)
- `-R`, `--dryrun`: run command without implementing any real actions
- `-S`, `--stats`: write csv file with summary statistics
- `-T`, `--threads`: number of threads to use in parallel mode
- `-v`, `--version`: show EXONtools current version
- `-W`, `--warnigns`: hide console output except warning and critical error messages

## Supplementary EXONtools commands

- `EXONtools.py force_install`: force PIP to install all required external python modules
- `EXONtools.py SetDefaultConfig`: write/overwrite default configuration settings in `dependencies.ini` file
- `EXONtools.py test`: run all unit tests

## Links to EXONtools dependencies

- [ABySS](https://github.com/bcgsc/abyss)
- [BBtools](https://jgi.doe.gov/data-and-tools/bbtools/)
- [BLAT](https://github.com/icebert/pblat)*
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [bwa](http://bio-bwa.sourceforge.net)
- [CAP3](http://seq.cs.iastate.edu/cap3.html)
- [CD-HIT](http://weizhongli-lab.org/cd-hit)
- [Cutadapt](https://cutadapt.readthedocs.io)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
- [FLASH](http://ccb.jhu.edu/software/FLASH)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [MultiQC](https://multiqc.info/)
- [SAMtools](http://www.htslib.org)
- [SPAdes](https://cab.spbu.ru/software/spades/)
- [TransABySS](https://github.com/bcgsc/transabyss)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Trinity](http://trinityrnaseq.github.io/)

\* rename `pblat` file to `blat` file to keep using the parallel BLAT

## Step-by-step pipeline example

```text
#! /usr/bin/env bash

##########################################
### Welcome to the EXONtools pipeline! ###
##########################################

### CREATE TEST DIRECTORY
mkdir TEST
cd TEST

#####################################
### STAGE A. Processing raw reads ###
#####################################

### STEP A1. Demultiplexing (demultiplex_reads)
../../src/EXONtools.py -S demultiplex_reads \
-R1 ../data/rnaseq/ST_rna_R1.fastq.gz \
-R2 ../data/rnaseq/ST_rna_R2.fastq.gz \
-o STEP_A1_DEMULTIPLEX -t 1  \
-b ../misc_files/rnaseq_indexes.dat --suffix "_demux"

### STEP A2. Formatting (format_reads)
../../src/EXONtools.py -T 2 -S format_reads \
-i STEP_A1_DEMULTIPLEX -o STEP_A2_PRECLEAN --rename --rqc

### STEP A3. Error correction (correct_reads)
../../src/EXONtools.py -T 2 -S correct_reads \
-i STEP_A2_PRECLEAN -o STEP_A3_CORRECT --suffix "_corr"



#####################################
### STAGE D. Processing raw reads ###
#####################################

### Concatenate gzip files to a single archive
cat ../data/exonseq/POPSET_S92_L004_R1.fastq.gz.a \
../data/exonseq/POPSET_S92_L004_R1.fastq.gz.b \
> POPSET_S92_L004_R1.fastq.gz
cat ../data/exonseq/POPSET_S92_L004_R2.fastq.gz.a \
../data/exonseq/POPSET_S92_L004_R2.fastq.gz.b \
> POPSET_S92_L004_R2.fastq.gz

### STEP D1. Demultiplexing (demultiplex_reads)
../../src/EXONtools.py -S demultiplex_reads \
-R1 POPSET_S92_L004_R1.fastq.gz -R2 POPSET_S92_L004_R2.fastq.gz \
-o STEP_D1_DEMULTIPLEX -t 1 \
-b ../misc_files/exon_indexes.dat --suffix "_demux"

### STEP D2. Formatting (format_reads)
../../src/EXONtools.py -T 2 -S format_reads \
-i STEP_D1_DEMULTIPLEX -o STEP_D2_PRECLEAN --rename --rqc

### STEP D3. Error correction (correct_reads)
### This should take some time for multiple libraries
### Use more threads to speed up
../../src/EXONtools.py -T 2 -S correct_reads \
-i STEP_A2_PRECLEAN -o STEP_D3_CORRECT --suffix "_corr"


```

## Future development

- speed up code implementation
- adding new EXONtools modules and/or dependencies
- replacing dependencies with EXONtools modules
- EXONtools docker image including all dependencies

## Glossary

- `ECS`: exon capture sequencing
- `FASTA`: text-based file format for storing nucleotide or peptide sequences
- `FASTQ`: text-based file format for storing nucleotide or peptide sequences produced by many next generation sequencing platforms
- `GBS`: genotyping by sequencing
- `Hybridization bait`: biotin-labeled oligonucleotide used to capture specific exon
- `NGS`: next generation sequencing
- `Non-model organism`: organism without high quality reference genome and genome annotation
- `PE reads`: paired-end reads
- `Pseudoreference`: annotated consensus assembly produced from various approaches/samples
- `SE reads`: single-end reads
- `SNP`: single nucleotide polymorphism
- `UCE`: ultraconserved element

## Citation

## Acknowledgements

## Contacts
