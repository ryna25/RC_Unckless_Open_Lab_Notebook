---
layout: post
title: PAML and preparing PAML output summary
date: 06 May 2024  
category: [ Computational Pipelines ]
tags: [ PAML, Comparative genomics ]
---

In this script, I'm doing PAML using the M0 model and Site model.

I'll generate my PAML control files for each ID using a for loop.


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

I'll run this file using the following two scripts - Script_PAML_1.sh & Script_PAML_2.sh 

Here is Script_PAML_1.sh 


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --partition=sixhour
#SBATCH --array=1-5000


echo "ID: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

module load paml
codeml "codeml-sites_$SLURM_ARRAY_TASK_ID".ctl
```

Here is Script_PAML_2.sh 


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --partition=sixhour
#SBATCH --array=5001-5759


echo "ID: $SLURM_ARRAY_TASK_ID" >> sh.$SLURM_JOBID.err

module load paml
codeml "codeml-sites_$SLURM_ARRAY_TASK_ID".ctl
```

Next, I'm going to store the omega (dN/dS) values from all PAML output files into a summary file using a for loop

I'm calling this script - GetPAMLsummary.sh


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
