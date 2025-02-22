# BCRA

Ball Covariance Ranking and Aggregation (BCRA) hypothesis test

Paper: Dai, W., & Zhang, H. (2025). Identifying genetic variants for brain connectivity using Ball Covariance Ranking and Aggregation. Journal of the American Statistical Association, 1–18. https://doi.org/10.1080/01621459.2025.2450837

## Code

Required codes are located under ./**code/\***, including `Ball-1.3.12.tar.gz`, `0-BCRA-transform-UKB.R` and `BallDistanceVector.cpp`.


## UKB Real Data Application Demonstration

### Data
*chr5_8_demo.RData*: a make-up data for showing how to apply subsample-BCRA on UKB as real genotype data contains sensitive information.

### Analysis Procedure

### Step 1: Calculate geodesic distance among functional connectivity matrices
Calculate geodesic distance (other distance measures can be applied depending on your interest) of functional connectivity matrices across each pair of subjects.

### Step 2: Partition Genotype into SNP-set
Partition Genotype into SNP-set: you can use gene/LD-block or physical locations (what we adopted) to do the partition. PLINK 2.0 software can be used to fulfill the partition.

### Step 3: Genotype QC on each SNP-set
Genotype QC: we exclude: 1) subjects with more than 10% missing genotypes; 2) variants with missing genotype rate larger than 10%; and 3) variants that failed the Hardy-Weinberg test at $10^{-6}$ level.

### Step 4: Run subsample-BCRA for each SNP-set
Run subsample-BCRA on a SNP-set using code located under `./demo/UKB/1-run-subsample-BCRA-SNPset.R`. Since we are unable to provide sensitive real data, we provided a demo dataset named as \*./data/chr5_8_demo.RData\* for illustration purpose to execute the code. The output of the code will give selected SNPs (`selected_snps_int`) and permutation p-value for this SNP-set `pval_results_perm`.

### Step 5: Summarize p-values associated with each SNP-set
We put a demonstration code `./demo/visualization/Figure_UKB_Pval.R` to plot p-values associated with each SNP-set (looking very similar to `./results/Figure_1.png`) with the data file `./demo/visualization/pval_UKB.RData` gives the SNPs selected with permuted p-value for each SNP-set presented in the manuscript.

### Step 6: Visualize the results
6. The connectivity plots in the manuscript (similar to `./results/Figure_2.png`) can be generated using `chordDiagram` function in R.

## Running Time and Resource Requirements on UKB

### a. Distance Calculation
- **Memory**: At least 100 GB per CPU
- **CPUs**: 4  
- **Time**: ~1 week  
- **Note**: Additional memory is required to write and store distance matrices.

### b. BCRA/Subsample Analysis
  - **Memory**: 32 GB per CPU 
  - **CPUs**: 4  
  - **Time**: BCRA: ~1 week per SNP-set with 5000 permutations. Subsample-BCRA: 2 days per SNP-set with 5000 permutations

=======
Ball Covariance Ranking and Aggregation (BCRA) hypothesis test
