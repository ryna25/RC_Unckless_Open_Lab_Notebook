---
layout: post
title: Preparing PAML input alignment and tree files
date: 06 May 2024  
category: [ Computational Pipelines ]
tags: [ PAML, Comparative genomics ]
---

The aim of this script is to get all the exon alignments and trees needed to do PAML over all genes on the X chromosome using *D. affinis* ST, SR1, SR2 and an outgroup *D. azteca*

I'm doing this on cluster in /work/unckless/a948g501/PAML/

Rob gave me a gtf file for all exon sequences for the *D. affinis* ST female genome. This file is called - Daffinis_STfemale_v5.1.ago2.gtf

This file has information for the whole genome transcriptome.

First I'm going to filter only te exons on ChrX and make a new gtf file for it.


```python
grep "ChrX_MullerAD" Daffinis_STfemale_v5.1.ago2.gtf > Daff_ST_ChrX.gtf
```

Rob gave me a *D. affinis* masked female genome file - Daffinis_STfemale_v5.1.masked.fasta

Now, I'm using gffread to create a multifasta file for my exons for all exons on ChrX.


```python
module load conda
conda activate gffread_env

gffread -w transcripts_X.fa -g Daffinis_STfemale_v5.1.masked.fasta Daff_ST_ChrX.gtf

conda deactivate
```

Next, I split my multi-fasta file into individual fasta files


```python
awk '/^>/{if(x>0) close(x".fasta");x++} {print > (x".fasta")}' transcripts_X.fa 
```

This created 5759 fasta files - each fasta file corresponding to one mRNA

Next, I obtain raw reads for *D. affinis* ST, SR1, and SR2 and *D. azteca*

I’ll copy D. affinis ST, SR1 and SR2 reads to this folder and concatenate them

I’m calling this script - GetDaffReads.sh


```python
#!/bin/bash

# Specify the directory containing your local files
File_dir="/work/unckless/a948g501/PAML/"

    #SR1
    scp a948g501@dtn.ku.edu:/resfs/GROUPS/MB/NIH0072614/Daffinis/RawReads/Fragments/Sample_54/54_GCCAAT_L001_R* ${File_dir}

    #SR2
    scp a948g501@dtn.ku.edu:/resfs/GROUPS/MB/NIH0072614/Daffinis/RawReads/Popgen/pool_57_*.fq.gz ${File_dir}

    #ST
    scp a948g501@dtn.ku.edu:/resfs/GROUPS/MB/NIH0072614/Daffinis/RawReads/Fragments/affinisF_frag_R*.fq.gz ${File_dir}

chmod 777 *.gz

#SR1 and SR2 reads have to be concatenated
cat 54_GCCAAT_L001_R1_001.fastq.gz 54_GCCAAT_L001_R1_002.fastq.gz > Daff_SR1_1.fastq.gz
cat 54_GCCAAT_L001_R2_001.fastq.gz 54_GCCAAT_L001_R2_002.fastq.gz > Daff_SR1_2.fastq.gz

gunzip -c Daff_SR1_1.fastq.gz > Daff_SR1_1.fastq
gunzip -c Daff_SR1_2.fastq.gz > Daff_SR1_2.fastq

cat pool_57_USPD16092225-N708-AK394_HV3F7BBXX_L2_1.fq.gz pool_57_USPD16092225-N708-AK394_HTYWGBBXX_L6_1.fq.gz > Daff_SR2_1.fastq.gz
cat pool_57_USPD16092225-N708-AK394_HV3F7BBXX_L2_2.fq.gz pool_57_USPD16092225-N708-AK394_HTYWGBBXX_L6_2.fq.gz > Daff_SR2_2.fastq.gz

gunzip -c Daff_SR2_1.fastq.gz > Daff_SR2_1.fastq
gunzip -c Daff_SR2_2.fastq.gz > Daff_SR2_2.fastq

#Remove SR1 and SR2 rawreads
rm 54_GCCAAT_L001_R*.fastq.gz
rm pool_57_*.fq.gz

#I'll change ST reads name to have similar nomenclature
gunzip -c affinisF_frag_R1.fq.gz > Daff_ST_1.fastq
gunzip -c affinisF_frag_R2.fq.gz > Daff_ST_2.fastq

rm Daff_SR*.fastq.gz
rm affinisF_frag_R*.fq.gz

```

Now I’ll save this script and make it executable and then run it using Script_GetDaffReads.sh


```python
chmod +x GetDaffReads.sh
```

Here is Script_GetDaffReads.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=100GB
#SBATCH --partition=sixhour
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


sh GetDaffReads.sh
```

I’ll copy the other reads from NCBI into this folder as well

I’m calling this script - GetOutgroupReads.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/PAML/

module load conda
conda activate sra-tool_env

#Copying D. azteca reads
fasterq-dump SRR12849542
mv SRR12849542_1.fastq Dazt_1.fastq
mv SRR12849542_2.fastq Dazt_2.fastq

conda deactivate

```

I’ll make this script executable and then run it using Script_GetOutgroupReads.sh


```python
chmod +x GetOutgroupReads.sh
```

Here is Script_GetOutgroupReads.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


sh GetOutgroupReads.sh
```

Now, we’re going to run fastqc on all our reads for a quality control check on our reads.

I’m calling this script - RunFastqc.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/PAML/

module load fastqc

file=$1

fastqc ${file}
```

I’ll make this script executable and run it using Script_RunFastqc.sh


```python
chmod +x RunFastqc.sh
```

Here is Script_RunFastqc.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=4GB
#SBATCH --partition=sixhour
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


sh RunFastqc.sh $1
```

To run this script, I used the following command to submit each fastqc job separately on cluster


```python
for file in *.fastq
do
sbatch Script_RunFastqc.sh $file
done
```

After looking at the results from fastqc, I’m trimming my reads using fastp and removing raw reads from this folder

I'm calling this script - RunFastp.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/PAML/

i=$1

module load conda
conda activate fastp

# Fastp to trim the reads
fastp -i "${i}_1.fastq" -I "${i}_2.fastq" -f 5 -t 5 -o "${i}_1.fastq.gz" -O "${i}_2.fastq.gz"
conda deactivate

# Remove raw reads
rm "${i}_1.fastq"
rm "${i}_2.fastq"
```

I’ll make this file executable and run it using Script_RunFastp.sh


```python
chmod +x RunFastp.sh
```

Here is Script_RunFastp.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --partition=sixhour
#SBATCH --array=1-
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1" >> sh.$SLURM_JOBID.err

sh RunFastp.sh $1
```

To run the Script, I want to use an iterative loop such that each i is submitted as a different job on cluster.

Here’s the command for it -


```python
for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dazt"
do 
sbatch Script_RunFastp.sh $i
done
```

Next, I'll align my files to each exon fasta file

I'm calling this script - BwaAlignIndex.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/PAML/

i=$1
j=$SLURM_ARRAY_TASK_ID


# Aligning to jth exon fasta file
module load bwa
module load samtools
bwa index "${j}".fasta
bwa mem "${j}".fasta "${i}_1.fastq.gz" "${i}_2.fastq.gz" | samtools view -hb -F 4 - | samtools sort - > "${i}_${j}.bam"

# Indexing the aligned file
samtools index "${i}_${j}.bam"
```

I'll make this script executable and run it using Script_BwaAlignIndex.sh


```python
chmod +x BwaAlignIndex.sh
```

Since I can only submit 5000 jobs at a time on cluster, I'll run it in 5 parts,

Here is Script_BwaAlignIndex_1.sh,
Script_BwaAlignIndex_2.sh,
Script_BwaAlignIndex_3.sh,
Script_BwaAlignIndex_4.sh,
Script_BwaAlignIndex_5.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --partition=sixhour
#SBATCH --array=1-1250
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1; j: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

sh BwaAlignIndex.sh $1
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --partition=sixhour
#SBATCH --array=1251-2500
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1; j: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

sh BwaAlignIndex.sh $1
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --partition=sixhour
#SBATCH --array=2501-3750
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1; j: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

sh BwaAlignIndex.sh $1
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --partition=sixhour
#SBATCH --array=3751-5000
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1; j: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

sh BwaAlignIndex.sh $1
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --partition=sixhour
#SBATCH --array=5001-5759
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1; j: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

sh BwaAlignIndex.sh $1
```

I'll run these jobs using for loop command


```python
for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dazt"
do 
sbatch Script_BwaAlignIndex_1.sh $i
done
```


```python
for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dazt"
do 
sbatch Script_BwaAlignIndex_2.sh $i
done
```


```python
for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dazt"
do 
sbatch Script_BwaAlignIndex_3.sh $i
done
```


```python
for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dazt"
do 
sbatch Script_BwaAlignIndex_4.sh $i
done
```


```python
for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dazt"
do 
sbatch Script_BwaAlignIndex_5.sh $i
done
```

Next I'll convert these bam files to fasta files to use for multiple sequence alignment.

I'll use these fasta files to generate a multiple sequence alignment using muscle and then make trees using iqtree.

I'm calling this script - bam2consensusfasta.sh


```python
#!/bin/bash

# Load modules 
module load samtools
module load bcftools
module load freebayes
module load conda
conda activate seqtk

j="$SLURM_ARRAY_TASK_ID"

for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dazt"
do 

    # Call SNPs on bam
    freebayes -f "${j}".fasta --ploidy 1 --bam "${i}_${j}".bam  > "${i}_${j}".vcf
    
    bgzip "${i}_${j}".vcf
    tabix "${i}_${j}".vcf.gz

    # Create consensus sequence from VCF
    bcftools consensus -f "${j}".fasta "${i}_${j}.vcf.gz" | awk -v id="${i}" 'BEGIN {FS = "\t"} {if (substr($1, 1, 1) == ">") {$1 = ">" id "_" substr($1, 2)} print}' > "${i}_${j}_consensus.fasta"

done

# Deactivate conda environment
conda deactivate

cat *"_${j}_consensus".fasta > "${j}_multi".fasta

#get the alignment
module load muscle
muscle -in "${j}_multi".fasta -out "alignment_${j}".phy -maxiters 1 -diags1

#make trees
module load iqtree
iqtree2 -s "alignment_${j}".phy -bb 1000 -wbt -nt AUTO
```

Now I'll make this script executable and run it using Script_bam2consensusfasta_1.sh and Script_bam2consensusfasta_2.sh


```python
chmod +x bam2consensusfasta.sh 
```

Here is Script_bam2consensusfasta_1.sh and Script_bam2consensusfasta_2.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --partition=sixhour
#SBATCH --array=1-5000


echo "j: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

sh bam2consensusfasta.sh
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --partition=sixhour
#SBATCH --array=5001-5759


echo "j: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

sh bam2consensusfasta.sh
```

Next, I'll do PAML using my alignment (.phy) and tree (.treefile) files.
