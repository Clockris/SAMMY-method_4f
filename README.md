# SAMMY-method 4f: bioinformatics analysis and plots generation

In this guide it is explained which are the bioinformatic steps to perform analysis of SAMMY-seq (4f protocol) as in "Biochemical properties of chromatin domains define genome compartmentalization" paper.

In this repository you are going to find four main folders: 
* `Data_processing`: where the scripts to perform the basic bioinformatics analysis starting from a fastq are stored
* `Input_samples`: where some input files for testing are stored (in this guide we are going to analyze them)
* `Analyzed`: where the results of the examples commands performed on `Input_samples` are stored 
* `Figure_making`: where the scripts to make the figures presented in the paper are stored.


## 1) Trimming

The first step of analysis is the trimming of the raw fastqs.
This operation is done through the Trimmomatic software and it has been implemented through the bash script `01-trimmer.sh`, which you can find in the `Data_processing` folder. For running this script the installation of Trimmomatic is required, we suggest to install through conda the version 0.39. 

The trimming script take as input:

1. the type of reads: it can be `single-end` or `paired-end`, according to the type of sequencing performed
2. the type of trimming to perform: it can be `SE` (single-end) or `PE` (paired-end), according to the type of input fastq provided
3. the number or cores that should be used for the analysis: an integer is required
4. the type of phred that will be used to asses the sequenced reads quality: it can be 'phred33' or 'phred64', depending on the Illumina pipeline used 
5. the name of the file where the log of Trimmomatic are to be saved
6. the name of the file where the statistics of Trimmomatic are to be saved
7. extra commands used by Trimmomatic (e.g. `'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15'`)
8. the path to the input file that should be trimmed
9. the name of the output trimmed file

You can run the script with a command as follows:
```
sh 01-trimmer.sh single-end SE 1 phred33 log_test stat_test 'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' input.fq.gz output.fq.gz
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
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20375_C2C12_2M_S2S_L003_R1_rand100000.fastq.gz \
	Analyzed/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.log \
	Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20376_C2C12_2M_S2L_L002_R1_rand100000.fastq.gz \
	Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimed.log \
	Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20377_C2C12_2M_S3_L003_R1_rand100000.fastq.gz \
	Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.log \
	Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/S20378_C2C12_2M_S4_L001_R1_rand100000.fastq.gz \
	Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz
```

## 2) Alignment
Once completed the trimming, the fastqs have to be aligned and then filtered to remove PCR duplicates, multimapping reads, bad quality aligned reads etc.
To perform this step we provide the script `02-aligner_and_filterer.sh`, which you can find in the `Data_processing` folder. To use this script, you need to install bwa (version 0.7.17-r1188), samtools (version 1.16.1) and picard (version 2.22) softwares. 

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
	Analyzed/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz Analyzed/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz Analyzed/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz Analyzed/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz \
	Analyzed \
	1540 \
	1
```


[^1] This index has been produced through the command `bwa index chr15.fa.gz` performed on chromosome chr15 of mm9 genome, dowloaded from UCSC at the following [link](https://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/chr15.fa.gz)
