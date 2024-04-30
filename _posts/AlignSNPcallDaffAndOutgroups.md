---
layout: post
title: Alignments of D. affinis and outgroup genome for SNP calling  
category: [ Computational Pipelines ]
tags: [ Sliding window trees, Comparative genomics ]
---

The aim of this pipeline is to call SNPs from *Drosophila affinis* SR1, SR2, ST and closely related species - *Drosophila algonquin, D. athabasca, D. azteca, and D. pseudoobscura*.

Rob gave me *D. affinis* SR1, SR2 and ST reads, and a reference masked genome. I am using short reads for other species from NCBI.

General steps in my pipeline after concatenating the reads - 

1. Run fastqc
2. Run fastp
3. Align using bwa
4. Index the bam file
5. Remove the reads, keep only aligned bam files and get stats on your aligned bam files
6. Call SNPs using bcftools-mpileup: separately for X(haploid) and autosomes(diploid)
7. Make trees across the X in sliding windows --> This is not a part of this pipeline - I'll do this in the next pipeline

I'm doing all of this on the cluster in /work/unckless/a948g501/SlidingTrees/

First I'm copying my masked genome file to this folder and indexing it - Daffinis_STfemale_v5.1.masked.fasta


```python
module load bwa 
bwa index Daffinis_STfemale_v5.1.masked.fasta

module load samtools
samtools faidx Daffinis_STfemale_v5.1.masked.fasta
```

Now I'll copy *D. affinis* ST, SR1 and SR2 reads to this folder and concatenate them

I'm calling this script - GetDaffReads.sh


```python
#!/bin/bash

# Specify the directory containing your local files
File_dir="/work/unckless/a948g501/SlidingTrees/"

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

Now I'll save this script and make it executable and then run it using Script_GetDaffReads.sh


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

I'll copy the other reads from NCBI into this folder as well 

I'm calling this script - GetOutgroupReads.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/SlidingTrees/

module load conda
conda activate sra-tool_env

#Copying D. algonquin Reads
fasterq-dump SRR5768634
mv SRR5768634_1.fastq Dalg_1.fastq
mv SRR5768634_2.fastq Dalg_2.fastq

#Copying D. athabasca reads
fasterq-dump SRR9997948
mv SRR9997948_1.fastq Datha_ea_1.fastq
mv SRR9997948_2.fastq Datha_ea_2.fastq

fasterq-dump SRR9997950
mv SRR9997950_1.fastq Datha_eb_1.fastq
mv SRR9997950_2.fastq Datha_eb_2.fastq

#Copying D. pseudoobscura reads
fasterq-dump SRR18151027
mv SRR18151027_1.fastq Dpse_1.fastq
mv SRR18151027_2.fastq Dpse_2.fastq

#Copying D. azteca reads
fasterq-dump SRR12849542
mv SRR12849542_1.fastq Dazt_1.fastq
mv SRR12849542_2.fastq Dazt_2.fastq

conda deactivate
```

I'll make this script executable and then run it using Script_GetOutgroupReads.sh


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

Now, we're going to run fastqc on all our reads for a quality control check on our reads.

I'm calling this script - RunFastqc.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/SlidingTrees/

module load fastqc

file=$1

fastqc ${file}
```

I'll make this script executable and run it using Script_RunFastqc.sh


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

After looking at the results from fastqc, I'm running the other steps in an iterative loop.

I'm calling this script - RunFastpBwaAlignIndex.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/SlidingTrees/

i=$1

module load conda
conda activate fastp

# Fastp to trim the reads
fastp -i "${i}_1.fastq" -I "${i}_2.fastq" -f 5 -t 5 -o "${i}_1.fastq.gz" -O "${i}_2.fastq.gz"
conda deactivate

# Aligning to the masked genome
module load bwa
module load samtools
bwa mem Daffinis_STfemale_v5.1.masked.fasta "${i}_1.fastq.gz" "${i}_2.fastq.gz" | samtools view -hb -F 4 - | samtools sort - > "${i}.bam"

# Indexing the aligned file
samtools index "${i}.bam"

# Stats on alignments
samtools flagstat "${i}.bam" > "${i}_bamStats.txt"
zcat "${i}"_*.fastq.gz | wc -l > "${i}_TotalLinesinReadsR1R2combined.txt"

# Remove all reads
rm "${i}_1.fastq"
rm "${i}_2.fastq"
rm "${i}_1.fastq.gz"
rm "${i}_2.fastq.gz"
```

I'll make this file executable and run it using Script_RunFastpBwaAlignIndex.sh


```python
chmod +x RunFastpBwaAlignIndex.sh
```

Here is Script_RunFastpBwaAlignIndex.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=100GB
#SBATCH --partition=sixhour
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err

echo "i: $1" >> sh.$SLURM_JOBID.err

sh RunFastpBwaAlignIndex.sh $1
```

To run the Script, I want to use an iterative loop such that each i is submitted as a different job on cluster.

Here's the command for it - 


```python
for i in "Daff_ST" "Daff_SR1" "Daff_SR2" "Dalg" "Datha_ea" "Datha_eb" "Dpse" "Dazt"
do 
sbatch Script_RunFastpBwaAlignIndex.sh $i
done
```

Next I want to combine the stats for my bam files into one file -


```python
for file in *.txt
do echo "==== $file ====" >> AllCombined_StatsBamFiles.txt
cat "$file" >> AllCombined_StatsBamFiles.txt
done
```

Now, I want to call SNPs from my aligned bam files

I'm using bcftools mpileup to call SNPs since these are divergent taxa!

Here's my Autosomes.bed file -


```python
Chr4_MullerB
Chr2_MullerE
Chr2.group4_MullerE
Chr3_MullerC
Chr5_MullerF
Unknown_69
mtDNA
```

I'm calling this script - CallSNPs_bcftools_mpileup.sh


```python
#!/bin/bash

cd /work/unckless/a948g501/SlidingTrees/

module load bcftools

# Haploid calling for ChrX_MullerAD
bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 5000 -f Daffinis_STfemale_v5.1.masked.fasta *.bam -r ChrX_MullerAD | bcftools call -vmO v --ploidy 1 -o SNPs_ChrX_haploid.vcf

# Diploid calling for autosomes
bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 5000 -f Daffinis_STfemale_v5.1.masked.fasta *.bam -R Autosomes.bed | bcftools call -vmO v -o SNPs_Autosomes_diploid.vcf

# SNP filtering and stats
bcftools filter -i 'QUAL > 30 && TYPE="snp"' SNPs_ChrX_haploid.vcf -o filtered_SNPs_ChrX_haploid.vcf
bcftools stats filtered_SNPs_ChrX_haploid.vcf > stats_filtered_SNPs_ChrX_haploid.txt

bcftools filter -i 'QUAL > 30 && TYPE="snp"' SNPs_Autosomes_diploid.vcf -o filtered_SNPs_Autosomes_diploid.vcf
bcftools stats filtered_SNPs_Autosomes_diploid.vcf > stats_filtered_SNPs_Autosomes_diploid.txt

```

I'll make this script executable and run it using Script_CallSNPs_bcftools_mpileup.sh


```python
chmod +x CallSNPs_bcftools_mpileup.sh
```

Here is Script_CallSNPs_bcftools_mpileup.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


sh CallSNPs_bcftools_mpileup.sh
```
