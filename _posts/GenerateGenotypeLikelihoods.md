I want to call SNPs from my alignments following mapping the reads from 90 male samples to the D. affinis female masked genome.

I am using angsd to call SNPs from these alignments.

I'm doing this on the cluster in /work/unckless/a948g501/PopGen/

I will first generate genotype likelihoods, and then call and filter SNPs using angsd. At the end, I will convert it into a vcf format using angsd.

I will do this separately for the X chromosome and then for the autosomes. I will also use this pipeline later for the Y chromosome.

I need to make a .bed file containing information about all other Chromosomes except the X.

I do this as follows - 


```python
module load samtools
samtools faidx Daffinis_STfemale_v5.1.masked.fasta 

cat Daffinis_STfemale_v5.1.masked.fasta.fai
```

I use the output from the cat .fai command to make the Autosomes.bed file


```python
nano Autosomes.bed
```


```python
Chr4_MullerB    1    49203969
Chr2_MullerE    1    44456471
Chr2.group4_MullerE    1    50000
Chr3_MullerC    1    23927763
Chr5_MullerF    1    1445284
Unknown_69    1    654907
mtDNA    1    15806
```

This is what goes into my Autosomes.bed file

Next, I generated a text file of my bam files list.

My pipeline to generate genotype likelihoods is called angsdGenerateGenotypeLikelihoods.sh and first I'm using this only for the X chromosome and then for the autosomes. Here is the script for angsdGenerateGenotypeLikelihoods.sh


```python
#!/bin/bash

module load angsd

ls /work/unckless/a948g501/PopGen/*.bam > bam_list.txt

# Specify the bed file containing regions to include
bed_file="Autosomes.bed"

# Generate genotype likelihoods
angsd -bam bam_list.txt -r ChrX_MullerAD -ref Daffinis_STfemale_v5.1.masked.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doHaploCall 2 -out output_prefix_X

angsd -bam bam_list.txt -rf $bed_file -ref Daffinis_STfemale_v5.1.masked.fasta -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 -out output_prefix_Autosomes

```

Next, I made my script executable using the following


```python
chmod +x angsdGenerateGenotypeLikelihoods.sh
```

Then, I ran the script on the cluster using the following Script_angsdGenerateGenotypeLikelihoods.sh


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH -o sh.%j.out
#SBATCH -e sh.%j.err


sh angsdGenerateGenotypeLikelihoods.sh
```


```python
sbatch Script_angsdGenerateGenotypeLikelihoods.sh
```
