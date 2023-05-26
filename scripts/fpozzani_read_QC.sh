#!/bin/bash

#### This script is intended to be the minimal command pipeline to execute all the assignment 
#### tasks from the beginning, as if nothing that follows had been done from command line 
#### before. 
#### This means this script is NOT intended to reproduce the assigment and should not be
#### executed unless all files, directories, and configurations are cleared.
#### Other than essential commands, I will only include commands that allow to answer
#### challenging questions in the assigment, and put my replies in comments.   



######################### Session 1 | Assignment setup

echo "script start: download and initial sequencing read quality control"
date

cd /shared/ifbstor1/projects/2314_medbioinfo/federico/MedBioinfo/analyses



######################### Session 2 | Programmatic access to NCBI SRA

### retrieve sequencing run ids
sqlite3 -batch -noheader -csv /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "select run_accession from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='fpozzani';" > fpozzani_run_accessions.txt

mkdir ../data/sra_fastq

ml sra-tools

### vdb configuaration
KEYS_DOWN=$'\e'[B
echo -n "C i T $KEYS_DOWN S o x" | vdb-config --interactive

### download FASTQ files via SRA suite API
cat fpozzani_run_accessions.txt | srun --cpus-per-task=1 --time=00:30:00 xargs fastq-dump --readids --gzip --outdir ../data/sra_fastq/ --disable-multithreading --split-e

### count the number of reads in a FASTQ file
zgrep -c '^@' ../data/sra_fastq/ERR6913192_1.fastq.gz
# reads are 322391


### how are the base call quality scores encoded?
zcat ../data/sra_fastq/ERR6913192_1.fastq.gz | head -n 20
# looks like Phred-33 ASCII encoding




######################### Session 3 | Manipulating raw sequencing FASTQ files with seqkit

ml seqkit/2.1.0

### use the appropriate seqkit sub-command to print statistics on each of the downloaded FASTQ files
srun --cpus-per-task=1 --time=00:05:00 seqkit stats -a --threads 1 ../data/sra_fastq/*

### compare the number of reads and bases to the values announced in the metadata provided by the authors
sqlite3 -batch /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "select * from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='fpozzani';"

### does-it look like the reads are un-trimmed (as produced straight out of the Illumina sequencer) or have been subjected to quality filtering/trimming ?
# it looks they were quality filtered since Q30 > 95% in all of them  


### can you use a seqkit sub-command to check if the FASTQ files have been de-replicated (duplicate identical reads removed) ?
### considering we will ultimately want to produce quantitative estimates of the pathogens present in the patient samples, which versions of the files is it better to work with ?
mkdir ../data/sra_fastq/duplicates
srun --cpus-per-task=3 --time=00:05:00 zcat ../data/sra_fastq/ERR6913310_2.fastq.gz | seqkit rmdup -s -o ../data/sra_fastq/duplicates/clean.fa.gz -d ../data/sra_fastq/duplicates/duplicated.fa.gz -D ../data/sra_fastq/duplicates/duplicated.detail.txt --threads 3
# in this file there were 78397 duplicated records. This means files have not been
# de-replicated. For a quantitative estimate, it's better to work with unprocessed
# non de-replicated files.  


###can you use a seqkit sub-command to guess if the FASTQ files have already been trimmed of their sequencing kit adapters
echo -e 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'>adaptors.txt
srun --cpus-per-task=3 --time=00:05:00 seqkit grep -s -f adaptors.txt ../data/sra_fastq/ERR6913310_2.fastq.gz -o adaptors_found.txt --threads 3

###start by looking for the whole adapters, then try shortening the adapters...
echo -e 'AGATCGGAAGAG\nAGATCGGAAGAG'>adaptors.txt
srun --cpus-per-task=3 --time=00:05:00 seqkit grep -s -f adaptors.txt ../data/sra_fastq/ERR6913310_2.fastq.gz -o adaptors_found.txt --threads 3
#looks like full lenght adaptors are not in sequences anymore, but shorter parts are.




######################### Session 4 | Quality control the raw sequencing FASTQ files with fastQC
mkdir fastqc
ml fastqc  

srun --cpus-per-task=2 --time=00:10:00 fastqc --outdir ./fastqc/ --threads 2 --noextract ../data/sra_fastq/ERR6913288_1.fastq.gz ../data/sra_fastq/ERR6913288_2.fastq.gz


srun --cpus-per-task=2 --time=00:30:00 xargs -I{} -a fpozzani_run_accessions.txt fastqc --outdir ./fastqc/ --threads 2 --noextract ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz




######################### Session 5 | Moving files from remote server to local laptop hard disk

### reads have been trimmed to exclude bases with low quality scores (often found nearer the end of reads)?
# yes
### reads have been trimmed to exclude sequencing library adapters?
# yes




######################### Session 7 | Merging paired end reads
ml flash2
mkdir ../data/merged_pairs

srun --cpus-per-task=2 flash2 --threads=2 -z --output-directory=../data/merged_pairs/ --output-prefix=ERR6913288.flash ../data/sra_fastq/ERR6913288_2.fastq.gz ../data/sra_fastq/ERR6913288_1.fastq.gz 2>&1 | tee -a fpozzani_flash2.log

### what proportion of your reads were merged successfully ?
# 84.81%
### use seqkit stat to check out the range of merged read lengths
srun --cpus-per-task=1 --time=00:05:00 seqkit stats -a --threads 1 ../data/merged_pairs/ERR6913288.flash.extendedFrags.fastq.gz

### check out the .histogram file : what does this suggest concerning the length of the DNA library insert sizes ?
vi ../data/merged_pairs/ERR6913288.flash.hist
# has info about how many y sequences are x bp long

srun --cpus-per-task=2 --time=00:30:00 xargs -a fpozzani_run_accessions.txt -n 1 -I{} flash2 --threads=2 -z --output-directory=../data/merged_pairs/ --output-prefix={}.flash ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a fpozzani_flash2.log

### compare how many base pairs you had in your initial unmerged reads, versus how many you have left after merging: have you lost information, or was it redundant information ?
srun --cpus-per-task=1 --time=00:05:00 seqkit stats -a --threads 1 ../data/merged_pairs/*.flash.extendedFrags.fastq.gz
srun --cpus-per-task=1 --time=00:05:00 seqkit stats -a --threads 1 ../data/merged_pairs/*.flash.notCombined*
# I have some information less, although I can't really say whether it was reduntant or not



######################### Session 8 | Use read mapping to check for PhiX contamination (and more...) 

mkdir ../data/reference_seqs

sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

export PATH=${HOME}/edirect:${PATH}
efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna

echo -e "export LC_CTYPE=en_US.UTF-8\nexport LC_ALL=en_US.UTF-8" >> ~/.bashrc
source ~/.bashrc

ml bowtie2 

mkdir ../data/bowtie2_DBs
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_DBs/PhiX_bowtie2_DB

mkdir analyses/bowtie

srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/PhiX_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S bowtie/fpozzani_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee bowtie/fpozzani_bowtie_merged2PhiX.log
###do you observe any hits against PhiX ?
#no, no reads aligned


efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fna

srun bowtie2-build -f ../data/reference_seqs/SC2_NC_045512.fna ../data/bowtie2_DBs/SC2_bowtie2_DB

srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/SC2_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S bowtie/fpozzani_merged2SC2.sam --threads 8 --no-unal 2>&1 | tee bowtie/fpozzani_bowtie_merged2SC2.log
### do you observe any samples which seem to have SC2 sequences ? 
# yes, 1313 reads aligned exactly 1 time. From .sam file, it looks like these come from 3 of my patients. 

######################### Session 9 | Combine quality control results into one unique report for all samples analysed

ml multiqc

srun multiqc --force --title "fpozzani sample sub-set" ../data/merged_pairs/ ./fastqc/ ./fpozzani_flash2.log ./bowtie/


date
echo "script end."
