---
layout: post
title: Align genomes for D. affinis SR1, SR2, ST
date: 13 April 2024    
category: [ Computational Pipelines ]
tags: [ Angsd SNP calling, Population genomics ]
---

This is Anjali's Notebook.

I am doing Population genomics analysis of D. affinis ST, SR1, SR2

Here is my script to cat reads, fastp, and align them to the genome

First, I indexed the fasta file of the masked genome on cluster


```python
module load bwa
bwa index Daffinis_STfemale_v5.1.masked.fasta
```

Rob gave me the female ST Daff genome assembly

This is the pipeline I used to cat reads, fastp them and align them to the genome on cluster.

The pipeline script is called DaffPopGenMappingPipeline.sh and is on my cluster in /work/unckless/a948g501/PopGen/

Below is the pipeline script: DaffPopGenMappingPipeline.sh 


```python
#!/bin/bash

# Specify the directory containing your local files
File_dir="/work/unckless/a948g501/PopGen/"
i=$1

    # Copy files from remote server to local directory
    scp a948g501@dtn.ku.edu:/resfs/GROUPS/MB/NIH0072614/Daffinis/RawReads/Popgen/pool_${i}_*.fq.gz ${File_dir} || { echo "Error: Failed to copy files for pool $i"; exit 1; }

    chmod 777 "${File_dir}pool_${i}_"*.fq.gz

    # Concatenate *1.fq.gz files
    cat "${File_dir}pool_${i}_"*1.fq.gz > "${File_dir}pool_${i}_1.fq.gz"

    # Concatenate *2.fq.gz files
    cat "${File_dir}pool_${i}_"*2.fq.gz > "${File_dir}pool_${i}_2.fq.gz"

    # Run fastp
    module load conda
    conda activate fastp
    fastp -f 5 -l 5 -i "${File_dir}pool_${i}_1.fq.gz" -I "${File_dir}pool_${i}_2.fq.gz" -o "${File_dir}fil.pool${i}_1.fq.gz" -O "${File_dir}fil.pool${i}_2.fq.gz"
    conda deactivate

    # Remove concatenated reads - only keep filtered reads
    rm "${File_dir}pool_${i}_"*.fq.gz 

    # Map to genome
    module load bwa
    module load samtools
    bwa mem Daffinis_STfemale_v5.1.masked.fasta "${File_dir}fil.pool${i}_1.fq.gz" "${File_dir}fil.pool${i}_2.fq.gz" | samtools view -hb -F 4 - | samtools sort - > "${File_dir}pool${i}.bam"
    samtools index "${File_dir}pool${i}.bam" 

    # Remove reads
    rm "${File_dir}fil.pool${i}_"*.fq.gz
```

After creating the script, I had to make it executable


```python
chmod +x DaffPopGenMappingPipeline.sh
```

To run this script using sbatch on cluster, the following sbatch script was used.

It is also in /work/unckless/a48g501/PopGen/

It is called Script_DaffPopGenMappingPipeline.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=100GB
#SBATCH --partition=sixhour
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1" >> sh.$SLURM_JOBID.err

sh DaffPopGenMappingPipeline.sh $1
```

To run the Script, I want to use an iterative loop such that each i is submitted as a different job on cluster.

I used the following command to execute my script.


```python
for i in $(seq 1 91)
do
sbatch Script_DaffPopGenMappingPipeline.sh $i
done
```

There is no file for i=73 or pool 73.

There are files for all other pools from i=1 to i=91
