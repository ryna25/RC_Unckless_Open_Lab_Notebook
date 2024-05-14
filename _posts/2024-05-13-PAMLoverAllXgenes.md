---
layout: post
title: PAML over all genes on X chromosome in D. affinis 
date: 13 May 2024  
category: [ Computational Pipelines ]
tags: [ PAML, Comparative genomics ]
---



The aim of this script is to get all the exon alignments and trees needed to do PAML over all genes on the X chromosome using D. affinis ST, SR1, SR2 and an outgroup D. azteca

I’m doing this on cluster in /work/unckless/a948g501/PAML/

Rob gave me a gtf file for all exon sequences for the D. affinis ST female genome. This file is called - Daffinis_STfemale_v5.1.ago2.gtf

This file has information for the whole genome transcriptome.


I will index my masked *D. affinis* female genome (Rob gave me the assembly)


```python
module load bwa
bwa index Daffinis_STfemale_v5.1.masked.fasta 
module load samtools
samtools faidx Daffinis_STfemale_v5.1.masked.fasta 
```

First I’m going to filter only the exons on ChrX and make a new gtf file for it.


```python
grep "ChrX_MullerAD" Daffinis_STfemale_v5.1.ago2.gtf > Daff_ST_ChrX.gtf
```

Now, I’m using gffread to create a multifasta file for my spliced exons for all exons on ChrX.


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

I'm copying alignments for *D. affinis* ST, SR1, SR2 and *D. algonquin*, *D. athabasca*, *D. azteca* and *D. pseudoobscura* from SlidingTrees folder to this folder.
Daff_SR1.bam      Daff_SR2.bam      Daff_ST.bam      Dalg.bam      Datha_ea.bam      Datha_eb.bam      Dazt.bam      Dpse.bam
Daff_SR1.bam.bai  Daff_SR2.bam.bai  Daff_ST.bam.bai  Dalg.bam.bai  Datha_ea.bam.bai  Datha_eb.bam.bai  Dazt.bam.bai  Dpse.bam.bai


For info on how I got these alignments done - look at https://anjaligupta1210.github.io/AG_Unckless_Open_Lab_Notebook/AlignSNPcallDaffAndOutgroups/

Next, I had to change permissions for all these bam files


```python
chmod 777 *.bam*
```

Next, I’ll align my files to each exon fasta file

I’m calling this script - AlignIndexForEachGene.sh


```python
#!/bin/sh
#
#SBATCH --job-name=auto.align.gene               
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=15               
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/PAML    
#SBATCH --mem-per-cpu=4gb            
#SBATCH --array=1-5759
#SBATCH --time=360

j=$SLURM_ARRAY_TASK_ID

module load bwa
module load samtools

bwa index "${j}".fasta

for i in "Daff_SR1" "Daff_SR2" "Daff_ST" "Dalg" "Datha_ea" "Datha_eb" "Dazt" "Dpse"
do
bwa mem "${j}".fasta "${i}".bam | samtools view -hb -F 4 - | samtools sort - > "${i}_${j}.bam"
samtools index "${i}_${j}.bam"
done

```

I'll make the script executable and then run it


```python
chmod +x AlignIndexForEachGene.sh
```

Next I’ll convert these bam files to consensus fasta files to use for multiple sequence alignment.

I’ll use these fasta files to generate a multiple sequence alignment using muscle and then make trees using iqtree.

I'm calling this script - bam2consensusfasta.sh


```python
#!/bin/bash
#
#SBATCH --job-name=auto.tree.alignment               
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=15               
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/PAML    
#SBATCH --mem-per-cpu=4gb            
#SBATCH --array=1-5759
#SBATCH --time=360


# Load modules 
module load samtools
module load bcftools
module load freebayes
module load conda
conda activate seqtk

j="$SLURM_ARRAY_TASK_ID"

for i in "Daff_SR1" "Daff_SR2" "Daff_ST" "Dalg" "Datha_ea" "Datha_eb" "Dazt" "Dpse"
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
iqtree2 -s "alignment_${j}".phy 
```

I'll make the script executable -


```python
chmod +x bam2consensusfasta.sh
```

I just wanted to check if I have alignments and treefiles for all transcripts and I do that here using this loop


```python
for ((j=1; j<=5759; j++)); do
    filename="alignment_${j}.phy.treefile"
    if [ ! -f "$filename" ]; then
        echo "Missing file: $filename"
        cat "alignment_${j}".phy
    fi
done
```

I’m doing PAML using the M0 model and Site model.

I’ll generate my PAML control files for each ID using a for loop.


```python
for i in {1..5759}
do
echo "      seqfile = "alignment_${i}".phy            * Path to the alignment file
     treefile = "alignment_${i}.phy".treefile           * Path to the tree file
      outfile = "out_sites_${i}".txt            * Path to the output file
   
        noisy = 3              * Display moderate amount of information on the screen
      verbose = 1              * Detailed output file

      seqtype = 1              * Codon data
        ndata = 1           * One gene alignment
        icode = 0              * Universal genetic code 
    cleandata = 0              * Do not remove sites with ambiguity data
		
        model = 0         * One ω for all branches (M0 and site models)
      NSsites = 0   1   2   7   8          * Heterogeneous ω Across Sites: Models M0 (0), M1a (1), M2a (2), M7 (7), and M8 (8)
    CodonFreq = 7        * Use mutation-selection model
      estFreq = 0        * Use observed frequencies to calculate fitness/freq pars
        clock = 0          * Assume no clock
    fix_omega = 0         * Enables option to estimate omega
        omega = 0.5        * Initial omega value" > codeml-sites_${i}.ctl
done
```

I’ll run these files for the PAML analysis


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --partition=sixhour
#SBATCH --array=1-5759


echo "ID: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

module load paml
codeml "codeml-sites_$SLURM_ARRAY_TASK_ID".ctl
```

Next, I’m going to store the omega (dN/dS) values from all PAML output files into a summary file using a for loop

I’m calling this script - GetPAMLsummary.sh


```python
#!/bin/bash

# Create a header for the output file
echo -e "File_ID\tGene_ID\tomega_dNdS" > PAML_output_summary.txt

for i in {1..5759}
do
    outfile="out_sites_${i}.txt"
    # Initialize variables to store grep results
    Gene_ID=$(grep ">Dazt_" "${i}_multi.fasta" | awk -F '>' '{sub("Dazt_", "", $2); print $2}')
    omega_dNdS=$(grep "omega (dN/dS) = " "$outfile" | awk '{print $NF}')

    # Append results to the output file
    echo -e "$i\t$Gene_ID\t$omega_dNdS" >> PAML_output_summary.txt

done
```

The output summary for PAML after running this script gets saved to PAML_output_summary.txt

The R file for the analysis is available at - https://github.com/anjaligupta1210/AG_Unckless_Open_Lab_Notebook/blob/master/rscripts/PAML.R
