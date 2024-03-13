# SAMMY-method 4f: bioinformatic analysis and plots generation

This repository contains the code (scripts) and a step by step tutorial explainign the bioinformatic steps to analyse SAMMY-seq data (4f protocol) as in the manuscript "Biochemical properties of chromatin domains define genome compartmentalization".

In this repository you are going to find four main folders: 
* `Scripts`: where the scripts to perform the basic bioinformatics analysis starting from FASTQ files are stored
* `Input_examples`: where some input files for testing are stored (in this guide we are going to analyze them)
* `Output_examples`: where the results of the examples commands performed on `Input_examples` are stored 

In this tutorial we are assuming that you will store all your outputs in a folder named `Output`. If you use as input the data stored in `Input_examples` double check your results with the ones stored in `Output_examples` to be sure everything worked fine.

## Before running the scripts
Before running the scripts described in this tutorial you should edit the `Scripts/environment_varialbles.src` with the path of the programs you are going to use and then run the following command (assuming you are located in the same folder where this README.md is placed):
```
source Scripts/environment_variables.src
```

## 1) Pre processing
### 1.1) Trimming

**Aim.** The first step of analysis is the trimming of the raw fastqs, this is done to remove from the reads all the adaptors and ensure a good standard quality for the alignment.

**Software.** This operation is done through the Trimmomatic software and it has been implemented through the bash script `1.1-trimmer.sh`, which you can find in the `Data_processing` folder. For running this script the installation of Trimmomatic is required, we suggest to install, through conda software, the version 0.39. 

**Input.** The trimming script take as input:
1. the type of reads: it can be `single-end` or `paired-end`, according to the type of sequencing performed
2. the type of trimming to perform: it can be `SE` (single-end) or `PE` (paired-end), according to the type of input fastq provided
3. the number or cores that should be used for the analysis: an integer is required
4. the type of phred that will be used to asses the sequenced reads quality: it can be 'phred33' or 'phred64', depending on the Illumina pipeline used 
5. the name of the file where the log of Trimmomatic are to be saved
6. the name of the file where the statistics of Trimmomatic are to be saved
7. extra commands used by Trimmomatic (e.g. `'ILLUMINACLIP:/path/to/TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15'`)
8. the path to the input file that should be trimmed
9. the name of the output trimmed file

**Output.** At the end of the step multiple output file will be created. One  containing the statistic of the trimming process (.stat), one where the log are stored (.log) and a fastq where the adapters from the original fastq have been removed plus a fastq containing the removed unpeared reads (only if two paired fastqs have been provided as input).

**Command.**

You can run the script with a command as follows:
```
sh 1.1-trimmer.sh single-end SE 1 phred33 log_test stat_test 'ILLUMINACLIP:/path/to/TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' input.fq.gz output.fq.gz
```

For the sake of this tutorial we are going to analyze the data in `Input_examples` through the following bash commands:
```
sh Scripts/1.1-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.log \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20375_C2C12_2M_S2S_L003_R1_rand100000.fastq.gz \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz

sh Scripts/1.1-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.log \
	Output/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20376_C2C12_2M_S2L_L002_R1_rand100000.fastq.gz \
	Input_examples/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz

sh Scripts/1.1-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.log \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20377_C2C12_2M_S3_L003_R1_rand100000.fastq.gz \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz

sh Scripts/1.1-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.log \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20378_C2C12_2M_S4_L001_R1_rand100000.fastq.gz \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz
```

### 1.2) Alignment
**Aim.** Once completed the trimming, the fastqs have to be aligned and then filtered to remove PCR duplicates, multimapping reads, bad quality aligned reads etc. 

**Software.** To perform this step we provide the script `1.2-aligner_and_filterer.sh`, which you can find in the `Scripts` folder. To use this script, you need to install bwa (version 0.7.17-r1188), samtools (version 1.17.1) and picard (version 2.22) softwares.

**Input.** The alignment and filtering script takes as input:
1. the number of cores that should be used for the analysis: an integer is required
2. the BWA indexed reference genome that will be used for the analysis
3. the path to the fastq file that should be aligned
4. the path to the folder where the output files will be stored
5. the numeric code that will be used by samtools -F command to filter reads (e.g 1540)
6. the minimum required quality for keeping reads during the filtering process (e.g. 1)

**Output.** From this analysis several bam files will be produced as output (coupled with the correspondent .bai indexes):
1. the file ending by "_full.bam" is the bam file just aligned to the reference genome
2. the file ending by "_mrkdup.bam" is the bam file aligned to the reference genome whose duplicates have been marked
3. the file ending by "_filtered.bam" is the bam file aligned to the reference genome whose reads have been filtered according to the input parameters

**Command.**
You can run the script with a command as follows:
```
sh 1.2-aligner_and_filterer.sh 4 mm9_UCSC_onlycanonical.fa.gz trimmed.fq.gz ./Output 1540 1
```

In this tutorial we are going to use the previously trimmed fastqs as input and the indexed genome `./Input_examples/Reference/chr15.fa.gz`[^1] as reference for this alignment process:
```
sh Scripts/1.2-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1

sh Scripts/1.2-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1

sh Scripts/1.2-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1

sh Scripts/1.2-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1
```

## 2) Post processing
The aligned and filtered reads produced by the previous commands are then processed to perform the analysis shown in the paper.

### 2.1) Coverage/single tracks
**Aim.** One of the most relevant information we could extract from the aligment and filtering of the sequencing files is how the reads distribute over the genome. This provide an information about how much the SAMMY-seq fraction in analysis is represented in each genomic region. The signal represented by this analysis is not the same that you can directly infer from the alignment (bam) file, but it is normalized by RPKM. This information is foundamental to call compartments.

**Software.** This analysis is performed throught the bamCoverage tool of DeepTools software suite (we suggest to use the version 3.4.3). To perform this step we provide the script `03a-genome_coverage__maker.sh` in `Scripts` folder.

**Input.** The script takes in input:
1. the bam file storing the information about the aligned and filtered reads created in the step 2) (e.g. input_filtered_aligned_trimmed.bam)
2. the output file where the reads will be stored (e.g. output.bw)
3. the output file type (e.g. bigwig)
4. the number of cores used to perform this process (e.g. 1)
5. the binsize for genome converage calculation expressed in base pairs (e.g. 5000)
6. the path to the sequences we want to exclude (e.g. blacklist.bed)
7. the genome effective size you can find on DeepTools *effective genome size* manual page (e.g. per mm9: 2620345972)
8. the method used for normalization (e.g. "RPKM")
9. the extension size of the aligned reads in base pairs (e.g. 250)

**Ouptut** The output is a bigwig file containing the reads distribution of each genomic region available in the input alignment file, normalized and without the blacklist regions.

**Command**
You can run the script as follow:
```
sh Scripts/2.1-genome_coverage__maker.sh input_filtered_aligned_trimmed.bam output.bw bigwig 1 5000 blacklist.bed 2620345972 RPKM 250
```
To test the script here provided you can run one or all of the several commands:

```
sh Scripts/2.1-genome_coverage__maker.sh \
	Output/02.3-01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S2S_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Scripts/2.1-genome_coverage__maker.sh \
	Output/02.3-01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S2L_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Scripts/2.1-genome_coverage__maker.sh \
	Output/02.3-01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S3_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Scripts/2.1-genome_coverage__maker.sh \
	Output/02.3-01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S4_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250
```

### 2.2) Comparisons/relative enrichement tracks
**Aim.** To increase the resolution of signal is often useful compare the signal coming from two different SAMMY-seq fractions, through this operation is possible better define regions of interests and check at the differences between two fractions. This operation is the same that is performed in ChIP-seq in performing a diffirence or a ratio between the signal coming from an IP and the one coming from an INPUT.

**Software.** To perform this operation we are using the Rscript `03b-comparison_maker.R` that you could find in the `Scripts` folder. For running this script `spp`, `data.table` and `Rcpp` libraries are required. We suggest to run the script using the R 3.5.2 version.

**Input.** For running this analysis/script four inputs are required:
1. the bam file of the first fraction we want compare (e.g. S2S.bam)
2. the bam file of the second fraction we want compare (e.g. S3.bam)
3. the chromosome size file, a file where for each lines there are three columns tab separated: (1) chromosome name as reported in the bams, (2) the coordinate of the first base of the chromosome (usually 1), (3) the coordinate of the last base of the chromosome. You can use as example the file `Input_examples/chr15.sizes`
4. the path and name.wig of the file where do you want to store your output (e.g. Output/S2SvsS3.wig)

**Ouput.** From running this script two (wiggle formatted) output files will results:
1. once ending with just ".wig" that is the output get.smoothed.tag.density(..., tag.shift = 0, background.density.scaling = TRUE) SPP function
2. once ending with just "_mle.wig" that is the output get.smoothed.enrichment.mle(..., tag.shift = 0, background.density.scaling = TRUE) SPP function

**Command.** 

You can run this script as follow:
```
Rscript 2.2-comparison_maker.R S2S.bam S3.bam chr.sizes S2SvsS3.wig
```

If you followed the previous command of this tutorial you can test the programm running the following line:
```
Rscript Scripts/2.2-comparison_maker.R Output/02.3-01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed_filtered.bam Output/02.3-01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed_filtered.bam Input_examples/chr15.size Output/03b-S2SvsS3_C2C12.wig
```

## 3) Downstream analyses and data visualization
### 3.0 ) Preparing file for compartment calling
**Aim.** This step is mandatory to split the output of 2.1 step (bigwig) per chromosome, in this was is it possible to call the compartments chromosome per chromosome.

**Software.** To perform this analysis the R software is required with installed the rtracklayer and data.table libraries, through this software is it possibile to run the `3.1-bw_binner.R` script (that leverage on the internal function sheet `3.1-bw_binner.Rl`)

**Input.** The input for this script are:
1. a bigwig generated from the 2.1 step (i.e Output_examples/03a-S4_C2C12.bw)
2. the information of the treatement of the sample (e.g. "untreated")
3. the information of the protocol performed (e.g. "4f")
4. the biological material on which the experiment has been perfomed (e.g. "C2C12")
5. the fraction in analysis (i.e. "S4")
6. the binning size: the size of the bin at which you want to rebin the sample (e.g. "250000")
7. the date in which the analysis has been performed (e.g. 20240213)
8. the reference genome name (e.g. "mm9")
9. the chromosome sizes of the chromosomes used formatted as a tsv where in the first field there is the chromosome name, in the second there is 1 and in the third there is the size of the chromosome (e.g. Input_examples/Split_chrs/mm9_chromsizes.bed)
10. the chromosome sizes of the chromosomes on which you want to call comparment formatted as 9 (it could be actually identical to/the same of 9). E.g. Input_examples/Split_chrs/mm9_only-canonical_chromsizes.bed.
11. the output folder where do you want to store the analysis (e.g. Output_examples/04-Splitted_chromosomes)

**Output.** The output will be a bigwig file for each chromosome presented in the 10th input element rebinned according to the 6th input parameter.

**Command.**
You can run this script as follow:
```
Rscript Scripts/3.1-bw_binner.R Output_examples/03a-S4_C2C12.bw untreated 4f C2C12 S4 250000 20240213 mm9 Input_examples/Split_chrs/mm9_chromsizes.bed Input_examples/Split_chrs/mm9_only-canonical_chromsizes.bed Output_examples/04-Splitted_chromosomes
```
If you followed the previous command of this tutorial you can test the program running the following line:
```
Rscript Scripts/3.1-bw_binner.R Output/03a-S2S_C2C12.bw untreated 4f C2C12 S2S 250000 20240213 mm9 Input_examples/Split_chrs/mm9_chromsizes.bed Input_examples/Split_chrs/mm9_only-canonical_chromsizes.bed Output_examples/04-Splitted_chromosomes

Rscript Scripts/3.1-bw_binner.R Output_examples/03a-S2L_C2C12.bw untreated 4f C2C12 S2L 250000 20240213 mm9 Input_examples/Split_chrs/mm9_chromsizes.bed Input_examples/Split_chrs/mm9_only-canonical_chromsizes.bed Output_examples/04-Splitted_chromosomes

Rscript Scripts/3.1-bw_binner.R Output_examples/03a-S3_C2C12.bw untreated 4f C2C12 S3 250000 20240213 mm9 Input_examples/Split_chrs/mm9_chromsizes.bed Input_examples/Split_chrs/mm9_only-canonical_chromsizes.bed Output_examples/04-Splitted_chromosomes

Rscript Scripts/3.1-bw_binner.R Output_examples/03a-S4_C2C12.bw untreated 4f C2C12 S4 250000 20240213 mm9 Input_examples/Split_chrs/mm9_chromsizes.bed Input_examples/Split_chrs/mm9_only-canonical_chromsizes.bed Output_examples/04-Splitted_chromosomes
```

### 3.1) Create compartment and subcompartment R objects

**Script.** `SAMMY_Compartment_getRobject_github.R`, `SAMMY_Compartment_utilityfunction_github.R`.

**Software.** `R`(suggested version 3.6.0) with the following packages: `data.table`, `GenomicRanges`, `pbapply`, `Matrix`.

**Aim.** Create and store all the required `.Rdata` objects used to call compartments by chromosome.
Namely, the following steps will be performed for each chromosome:
1. Read and save 1-dimensional (1D) tracks (i.e., S2S, S2L, S3 and S4 fractions);
2. Read and save 2-dimensional (2D) data (i.e., Hi-C matrix);
3. Read 1D and 2D data, compute and save (a) correlation matrix and (b) bin coordinates;
4. Read correlation matrices, compute eigenvectors and store the first principal component (PC1);
5. *(optional)* Read bin coordinates and compute the chromHMM states occupancy for each bin.

> [!NOTE]
> For the compartments estimation in [*"Compute A and B compartments"*](#32-compute-a-and-b-compartments) all the steps (1, 2, 3a, 3b and 4) and the created objects are required. For subcompartment estimation in [*"Compute subcompartments"*](#33-compute-subcompartments) only step 1 is required.



**Input files.**
- SAMMY-seq bigwig files, split by chromosome (see output in [*"Preparing file for compartment calling"*](#30--preparing-file-for-compartment-calling));
- bins table and Hi-C cooler files split by chromosome (recommend resolution 250 kb);
- *(optional)* ChormHMM segmentation of the genome (.dense file);

> [!TIP]
> Bins table and Hi-C cooler file can be obtained with the following command ([cooler](https://cooler.readthedocs.io/en/latest/datamodel.html)):
> ```
> cooler dump --table bins --header --out HiC_250000_bins.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/250000
> cooler dump --range chr18:0-80373285 --balanced --header --out HiC_250000_chr18_counts.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/250000
> ```
> Example mcool file can be downloded from [4DN portal](https://data.4dnucleome.org/files-processed/4DNFIMDOXUT8/).


**Input variable.** Input variables can be set within the script as follow:
```
# <chr> Path to the main folder
prefix.path<-"/Analysis/" 
# <chr> Name used as prefix for the created RData objects
prefix.name<-"Fibroblast_paper_202301_"

# <chr> Bin size in bp used to create the bigwig input files
binsize<-c("250000")
# <chr> Chromosome to be analyzed
CHR<-c("chr18")

# <chr> Path to the folder that contains the Hi-C cooler files
dir.hic<-paste0(prefix.path,"Data/4DN/")
#<chr> (optional) Path to the folder that contains ChromHMM dense file 
dir.hmm<-paste0(prefix.path,"Data/ChromHMM/")
```

**Output.** The created .RData objects will be stored in the main folder Robject/ and the correspondent subfolders. Namely, each chromosome file will contain:
1. `Track/`: list of GenomicRanges objects.
2. `HiC/`: list of matrix.
3. a)`Correlation/`: list of Matrix; and b) `Bin/`: GenomicRanges object.
4. `Eigenvector/`: list of vector.
5. *(optional)* `Bin/`: GenomicRanges object.


### 3.2) Compute A and B compartments

**Script.** `SAMMY_Compartment_Analysis_github.R`, `SAMMY_Compartment_utilityfunction_github.R`.

**Software.** `R`(suggested version 3.6.0) with the following packages: `data.table`, `GenomicRanges`, `Matrix`.

**Aim.** Get A and B compartments for SAMMY-seq and Hi-C samples, compare them and summarize the results.

**Input files.**
- `Bin/` R objects (see output 3b, in [*"Create R object"*](#31-create-compartment-and-subcompartment-r-objects));
- `Eigenvector/` R objects (see output 4, in [*"Creare R object"*](#31-create-compartment-and-subcompartment-r-objects)).

**Input variables.** Input variables can be set within the script as follow:
```
# <chr> Bin size in bp used to create the bigwig input files and the Hi-C cooler files
binsize<-"250000"
# <chr> Chromosome to be analyzed
CHR<-"chr18"
# <boolean> If true all eigenvector will be divide by the abs(max) (visual purpose) 
scale.egv<-TRUE
# <chr> Name of the data used as reference, to compare with
ref.hic<-"HiC"
```

**Output.** A bed file with genomic coordinates for each bin with a color indication of A (`#ffb344`) and B (`#00a19d`) compartment, for each experiment/sample. The output files will be stored in the subfolder `Results/Compartment/`. Summarizing plots will be saved in the folder `Plot/`.



### 3.3) Compute subcompartments

**Script.** `SAMMY_Subcompartment_Calder_Analysis_github.R`, `SAMMY_Subcompartment_Calder_UtilityFunctions_github.R`.

**Software.** `R` (suggested version 3.6.0) with the following packages: `data.table`, `GenomicRanges`, `doParallel`, `R.utils`, `factoextra`, `maptools`, `CALDER` (version 1.0, 2020-09-01), `patchwork`, `dendextend`, `scales`.
> [!NOTE]
> In this step, we re-implemented some of the CALDER algorithms to accommodate the SAMMY-seq data format. We still keep the core set of functions (`remove_blank_cols`, `fast_cor`, `generate_compartments_bed`, `HighResolution2Low_k_rectangle`, `get_PCs`, `bisecting_kmeans`, `project_to_major_axis`, `get_best_reorder`, `get_cluser_levels`), with the default parameters and the steps intended in the [original method](https://github.com/CSOgroup/CALDER).


**Aim.** Convert Hi-C cooler file in CALDER format and get the eigth sub-compartments segmentation as described in CALDER, for each SAMMY-seq and Hi-C sample.

**Input files.**
- SAMMY-seq bigwig files split for each chromosome (see point 3.1. [*"Create R object"*](#31-create-compartment-and-subcompartment-r-objects));
- bins table and Hi-C cooler files split by chromosome (recommend resolution 50 kb).


> [!TIP]
> Bins table and Hi-C cooler file can be obtained with the following command ([cooler](https://cooler.readthedocs.io/en/latest/datamodel.html)):
> ```
> cooler dump --table bins --header --out HiC_50000_bins.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/50000
> cooler dump --range chr18:0-80373285 --balanced --header --out HiC_50000_chr18_counts.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/50000
> ```
> Example mcool file can be downloded from [4DN portal](https://data.4dnucleome.org/files-processed/4DNFIMDOXUT8/).


**Input variables.** Input variables can be set within the script as follow:
```
# <chr> Path to the main folder
prefix.path<-"/Analysis/" 
# <chr> Name used as prefix for the created RData objects
prefix.name<-"Fibroblast_paper_202301_"

# <chr> Track directory (see point 3.0)
dir.track<- paste0(prefix.path,"Robjects/Track/")
# <chr> Path to the folder that contains the Hi-C cooler files
dir.hic<-paste0(prefix.path,"Data/4DN/")

# <chr> Bin size in bp used to create the bigwig input files and the Hi-C cooler files
binsize<-"50000"
```


**Output.** A bed file with genomic coordinates for each bin with color indication of A.1.1 (most open, `#1d4f60`), A.1.1. (`#1d4f60`), A.1.2 (`#36877a`), A.2.1 (`#6dbc90`), A.2.2 (`#c4e6c3`), B.1.1 (`#ecda9a`), B.1.2 (`#f3ad6a`) ,B.2.1 (`#f97b57`) and B.2.2. (`#ee4d5a`, most closed) subcompartment, for each experiment/sample. The output files will be stored in the sub folder `Results/Subcompartment/`. Summarizing plots will be saved in the folder `Plot/`.




[^1]: This index has been produced through the command `bwa index chr15.fa.gz` performed on chromosome chr15 of mm9 genome, dowloaded from UCSC at the following [link](https://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/chr15.fa.gz).
