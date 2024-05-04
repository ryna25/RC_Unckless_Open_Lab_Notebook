---
layout: post
title: Make trees in sliding windows
date: 04 May 2024  
category: [ Computational Pipelines ]
tags: [ Sliding window trees, Comparative genomics ]
---

Navigate to /work/unckless/a948g501/SlidingTrees/ on the cluster

I'm using Devon's pipelines to make trees from my SNP vcf files.

I need an all sites vcf file rather than just SNPs so, I'm going to filter my raw vcf file for QC but use all called sites including invariants included


```python
bcftools filter -i 'QUAL > 30' SNPs_ChrX_haploid.vcf -o filtered_allsites_ChrX_haploid.vcf

bcftools filter -i 'QUAL > 30' SNPs_Autosomes_diploid.vcf -o filtered_allsites_Autosomes_diploid.vcf
```

The first step is to take the vcf files and convert it into a phylip file so that it can just be treated like any other DNA sequence alignment (although a rather large one).

Here's the code to do this -


```python
##get conversion script
wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py

##load python
module load python

##convert file
python vcf2phylip.py -i /work/unckless/a948g501/SlidingTrees/filtered_allsites_ChrX_haploid.vcf
python vcf2phylip.py -i /work/unckless/a948g501/SlidingTrees/filtered_allsites_Autosomes_diploid.vcf
```

The PHYLIP matrix is saved to: filtered_allsites_ChrX_haploid.min4.phy and filtered_allsites_Autosomes_diploid.min4.phy

Now I will cut up the phylip

Credit to Jonothan Chang for a very clear explanation and framework for how best to do this (see: https://jonathanchang.org/blog/splitting-a-concatenated-raxml-style-phylip-file/)


For filtered_allsites_ChrX_haploid.min4.phy -

Length of the phylip alignment is 3,069,347

So I will generate 306 10Kb alignments and ignore the last 9.347Kb

To do so, I will define two variables in a for loop.

One (i) will track 1-306, over the 306 iterations
one (j) will track 1-3,050,001 in increments of 10Kb, over the 306 iterations

The code '$((j+9999))' means add 9,999 to j each iteration, which will track 10,000-3,060,000 in increments of 10Kb, over the 306 iterations

Run this code to generate the necessary partition file


```python
j=1
for i in {1..306}
do
echo "DNA, p${i}=${j}-$((j+9999))" >> partitions_X.txt
j=$((j+10000))
done
```

For filtered_allsites_Autosomes_diploid.min4.phy -

Length of the phylip alignment is 6,033,053

So I will generate 603 10Kb alignments and ignore the last 3.053Kb

To do so, I will define two variables in a for loop.

One (i) will track 1-603, over the 603 iterations
one (j) will track 1-6,020,001 in increments of 10Kb, over the 603 iterations

The code '$((j+9999))' means add 9,999 to j each iteration, which will track 10,000-6,030,000 in increments of 10Kb, over the 603 iterations

Run this code to generate the necessary partition file


```python
j=1
for i in {1..603}
do
echo "DNA, p${i}=${j}-$((j+9999))" >> partitions_Autosome.txt
j=$((j+10000))
done
```

We will download the necessary python script to unconcatenate the phylip based on the above partition file


```python
curl -LO https://gist.githubusercontent.com/jonchang/34c2e8e473ec2e8f50574671e62c3365/raw/unconcatenate_phylip.py
```

We will run the python script on this file to cut it up and read each subset file out into a new directory


```python
mkdir partitioned_X.phylips
module load python
python unconcatenate_phylip.py filtered_allsites_ChrX_haploid.min4.phy partitions_X.txt --prefix=partitioned_X.phylips/
```


```python
mkdir partitioned_Autosome.phylips
module load python
python unconcatenate_phylip.py filtered_allsites_Autosomes_diploid.min4.phy partitions_Autosome.txt --prefix=partitioned_Autosome.phylips/
```

The next step is just making a tree for each phylip alignment.

We are going to do this by submitting a 'batch job' on the cluster that spawns a unique job for each phylip alignment.

Here is my script for Script_RunGeneTreeArray_X.sh


```python
#!/bin/sh
#
#SBATCH --job-name=auto.gene.trees.X               
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=15               
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/SlidingTrees/partitioned_X.phylips    
#SBATCH --mem-per-cpu=2gb            
#SBATCH --array=1-306
#SBATCH --time=360


#To be able to later use bootstrapping in ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option -B 
#we do not specify a substition model with the -m option and therefore allow IQ-TREE to automatically select the best-fitting model.
#we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option --wbt.
#-nt AUTO ensures that all available CPU's are used

##load module iqtree
module load iqtree
iqtree2 -s /work/unckless/a948g501/SlidingTrees/partitioned_X.phylips/DNA_p$SLURM_ARRAY_TASK_ID.phylip -bb 1000 -wbt -nt AUTO
```

Make the Script executable -


```python
chmod +x Script_RunGeneTreeArray_X.sh
```

Here is my script for Script_RunGeneTreeArray_Autosome.sh


```python
#!/bin/sh
#
#SBATCH --job-name=auto.gene.trees.Autosome              
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=15              
#SBATCH --partition=sixhour            
#SBATCH --chdir=/work/unckless/a948g501/SlidingTrees/partitioned_Autosome.phylips    
#SBATCH --mem-per-cpu=2gb           
#SBATCH --array=1-603
#SBATCH --time=360


#To be able to later use bootstrapping in ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option -B 
#we do not specify a substition model with the -m option and therefore allow IQ-TREE to automatically select the best-fitting model.
#we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option --wbt.
#-nt AUTO ensures that all available CPU's are used

##load module iqtree
module load iqtree
iqtree2 -s /work/unckless/a948g501/SlidingTrees/partitioned_Autosome.phylips/DNA_p$SLURM_ARRAY_TASK_ID.phylip -bb 1000 -wbt -nt AUTO
```

Make the Script executable -


```python
chmod +x Script_RunGeneTreeArray_Autosome.sh
```

This will create a lot of output files but the only thing I really need to worry about is the output tree file. 
Once I have all of the output trees, I use the cat command to paste them all into a single file. 


```python
cat partitioned_X.phylips/DNA_p*.phylip.treefile > OutputTrees_X.trees
```


```python
cat partitioned_Autosome.phylips/DNA_p*.phylip.treefile > OutputTrees_Autosome.trees
```

Now we can calculate the proportion of trees supporting various topologies and visualize them using R 
