#!/bin/bash

# change directory
# This is where the raw fast files are. All of them are in the same folder 
cd /tgen_labs/jfryer/projects/chemobrain/1MPI/bulkRNA/

# create file with list of R1 samples. Data is paired end R1 and R2. 
# We only need to collect the read information once per sample. The read information is in both the R1 and R2 fastq files. 
ls -1 | grep _R1_ > R1Samples_1MPI.txt

# loops through list 
touch sampleReadInfo_1MPI.txt # creates an empty file
for sample in `cat R1Samples_1MPI.txt`; do
    zcat ${sample} | head -1 >> sampleReadInfo_1MPI.txt # read the first line of each R1 file
done;

# mv the files 
mv R1Samples_1MPI.txt  /tgen_labs/jfryer/kolney/chemobrain/scripts/01_processing/R1Samples_1MPI.txt
mv sampleReadInfo_1MPI.txt /tgen_labs/jfryer/kolney/chemobrain/scripts/01_processing/sampleReadInfo_1MPI.txt

cd /tgen_labs/jfryer/kolney/chemobrain/scripts/01_processing/
paste -d "\t" R1Samples_1MPI.txt sampleReadInfo_1MPI.txt > sampleReadGroupInfo_1MPI.txt # create sample info
rm R1Samples_1MPI.txt
rm sampleReadInfo_1MPI.txt