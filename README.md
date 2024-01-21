# SAMMY-method 4f: bioinformatics analysis and plots generation

This tutorial illustrate the use of the bioinformatic pipeline to analyze SAMMY-seq (4f protocol) data as presented in the "Biochemical properties of chromatin domains define genome compartmentalization" manuscript. The code in this repository also contains scripts to obtain the figures of the same manuscript.

In this repository you are going to find four main folders: 
* `Data_processing`: where the scripts to perform the basic bioinformatics analysis starting from a fastq are stored
* `Input_samples`: where some input files for testing are stored (in this guide we are going to analyze them)
* `Analyzed`: where the results of the examples commands performed on `Input_samples` are stored 
* `Figure_making`: where the scripts to make the figures presented in the paper are stored.

## Before run the scripts
Before run the script described in this tutorial you should edit the `Data_processing/environment_varialbles.src` with the path of the programs you are going to use and then run the following command (assuming you are located in the same folder where this README.md is placed):
```
source Data_processing/environment_variables.src
```

## 1) Trimming

The first step of analysis is the trimming of the raw fastqs.
This operation is done through the Trimmomatic software and it has been implemented through the bash script `01-trimmer.sh`, which you can find in the `Data_processing` folder. For running this script the installation of Trimmomatic is required, we suggest to install, through conda software, the version 0.39. 

The trimming script take as input:

1. the type of reads: it can be `single-end` or `paired-end`, according to the type of sequencing performed
2. the type of trimming to perform: it can be `SE` (single-end) or `PE` (paired-end), according to the type of input fastq provided
3. the number or cores that should be used for the analysis: an integer is required
4. the type of phred that will be used to asses the sequenced reads quality: it can be 'phred33' or 'phred64', depending on the Illumina pipeline used 
5. the name of the file where the log of Trimmomatic are to be saved
6. the name of the file where the statistics of Trimmomatic are to be saved
7. extra commands used by Trimmomatic (e.g. `'ILLUMINACLIP:/path/to/TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15'`)
8. the path to the input file that should be trimmed
9. the name of the output trimmed file

You can run the script with a command as follows:
```
sh 01-trimmer.sh single-end SE 1 phred33 log_test stat_test 'ILLUMINACLIP:/path/to/TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' input.fq.gz output.fq.gz
```

For the sake of this tutorial we are going to analyze the data in `Input_samples` through the following bash commands:
```
sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.log \
	Analyzed/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20375_C2C12_2M_S2S_L003_R1_rand100000.fastq.gz \
	Analyzed/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.log \
	Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20376_C2C12_2M_S2L_L002_R1_rand100000.fastq.gz \
	Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.log \
	Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20377_C2C12_2M_S3_L003_R1_rand100000.fastq.gz \
	Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.log \
	Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20378_C2C12_2M_S4_L001_R1_rand100000.fastq.gz \
	Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz
```

## 2) Alignment
Once completed the trimming, the fastqs have to be aligned and then filtered to remove PCR duplicates, multimapping reads, bad quality aligned reads etc.
To perform this step we provide the script `02-aligner_and_filterer.sh`, which you can find in the `Data_processing` folder. To use this script, you need to install bwa (version 0.7.17-r1188), samtools (version 1.17.1) and picard (version 2.22) softwares. 

The alignment and filtering script takes as input:

1. the number of cores that should be used for the analysis: an integer is required
2. the BWA indexed reference genome that will be used for the analysis
3. the path to the fastq file that should be aligned
4. the path to the folder where the output files will be stored
5. the numeric code that will be used by samtools -F command to filter reads (e.g 1540)
6. the minimum required quality for keeping reads during the filtering process (e.g. 1)

You can run the script with a command as follows:
```
sh 02-aligner_and_filterer.sh 4 mm9_UCSC_onlycanonical.fa.gz trimmed.fq.gz ./Analyzed 1540 1
```

In this tutorial we are going to use the previously trimmed fastqs as input and the indexed genome `./Input_sample/Reference/chr15.fa.gz`[^1] as reference for this alignment process:
```
sh Data_processing/02-aligner_and_filterer.sh \
	1 \
	Input_samples/Reference/chr15.fa.gz \
	Analyzed/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz \
	Analyzed \
	1540 \
	1

sh Data_processing/02-aligner_and_filterer.sh \
	1 \
	Input_samples/Reference/chr15.fa.gz \
	Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz \
	Analyzed \
	1540 \
	1

sh Data_processing/02-aligner_and_filterer.sh \
	1 \
	Input_samples/Reference/chr15.fa.gz \
	Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz \
	Analyzed \
	1540 \
	1

sh Data_processing/02-aligner_and_filterer.sh \
	1 \
	Input_samples/Reference/chr15.fa.gz \
	Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz \
	Analyzed \
	1540 \
	1
```
## 3) Post processing
The aligned and filtered reads produced by the previous commands are then processed to perform the analysis shown in the paper.

### 3.1 Coverage/single tracks
From the alignment files the reads distribution along the reference genome is calculated throught the bamCoverage tool of DeepTools software suite and the information is stored in a file. To perform this step we provide the script `03a-genome_coverage__maker.sh` in `Data_processing` folder.

The script takes in input:
1 - the bam file storing the information about the aligned and filtered reads created in the step 2) (e.g. input_filtered_aligned_trimmed.bam)
2 - the output file where the reads will be stored (e.g. output.bw)
3 - the output file type (e.g. bigwig)
4 - the number of cores used to perform this process (e.g. 1)
5 - the binsize for genome converage calculation expressed in base pairs (e.g. 5000)
6 - the path to the sequences we want to exclude (e.g. blacklist.bed)
7 - the genome effective size you can find on DeepTools effective genome size manual page (e.g. per mm9: 2620345972)
8 - the method used for normalization (e.g. "RPKM")
9 - the extension size of the aligned reads in base pairs (e.g. 250)

You can run the script as follow:
```
sh Data_processing/03a-genome_coverage__maker.sh input_filtered_aligned_trimmed.bam output.bw bigwig 1 5000 blacklist.bed 2620345972 RPKM 250
```
To test the script here provided you can run one or all of the several commands:

```
sh Data_processing/03a-genome_coverage__maker.sh \
	Analyzed/02.3-01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed_filtered.bam \
	Analyzed/03a-S2S_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_samples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Data_processing/03a-genome_coverage__maker.sh \
	Analyzed/02.3-01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed_filtered.bam \
	Analyzed/03a-S2L_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_samples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Data_processing/03a-genome_coverage__maker.sh \
	Analyzed/02.3-01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed_filtered.bam \
	Analyzed/03a-S3_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_samples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Data_processing/03a-genome_coverage__maker.sh \
	Analyzed/02.3-01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed_filtered.bam \
	Analyzed/03a-S4_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_samples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250
```

[^1]: This index has been produced through the command `bwa index chr15.fa.gz` performed on chromosome chr15 of mm9 genome, dowloaded from UCSC at the following [link](https://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/chr15.fa.gz)
