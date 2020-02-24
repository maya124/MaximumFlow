# MaximumFlow

## Maximum Flow Formulation Identifies High-Confidence Variants in Simple Repeat Sequences Associated with Autism Spectrum Disorder
#### Maya Varma

This repository contains the source code for implementing maximum flow for dimensionality reduction and improved feature stability in high-dimensional whole genome datasets. The algorithm uses the presence of linkage disequilibrium (LD) to cluster similar variants and reduce the size of the feature space.

## Implementation

This code base is implemented in Python 3.6.1 and requires the following libraries: numpy (1.16.4). Other necessary libraries include: bcftools (1.8), plink2 (2.0).


Begin by performing k-fold cross-validation with any l1-regularized machine learning algorithm on a high-dimensional whole genome dataset. In this implementation, we use k=5. Isolate all predictive variants (nonzero coefficient scores) from each fold and create an input folder with k files, each in the following tab-delimited format:  

```
chr-variant1	coefficient_score
chr-variant2	coefficient_score
...
```

First, we will generate a formatted list of all variants from the cross fold analysis. Run the following script:
```
python step_0.py  <input_folder>
```

This will generate two files in the input folder: variant_ranges.tsv and vcf_input.txt. Next, run the following command to filter these variants from whole genome vcf files:
```
bcftools view ${vcf_filename} -o ${outputfilename} --output-type z --regions-file vcf_input.txt
```

Now, we use plink to generate the necessary files for computing LD. This should create a series of files, include .bim and .bed files. 
```
plink2 --vcf $vcf_filename --extract range variant_ranges.tsv --make-bed --out $prefix
```

Run the following script to assign id numbers to variants so that LD scores can be computed. Specify the prefix of the .bim and .bed files (defined in the previous command).
```
python step_1.py <prefix>
```

Next, we use this information to calculate pairwise LD scores for variants across different folds. $num1 and $num2 represent the fold numbers to be compared (e.g. 1 and 3). 
```
srun ipython step_2.py $input_folder $prefix $num1 $num2
```

Finally, run the maximum flow algorithm to identify clustered groups of variants in LD. 
```
python maxflow.py 1_2_3_4_5
```

Results are stored in "results.txt" and "path.txt".


## Acknowledgments

* Maximum flow set-up and Ford-Fulkerson code adapted from: https://www.geeksforgeeks.org/ford-fulkerson-algorithm-for-maximum-flow-problem/

