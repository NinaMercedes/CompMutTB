# CompMut-TB - A framework to identify compensatory mutations from *Mycobacterium tuberculosis* whole genome sequences
## Setup
To setup and install dependencies for CompMutTB, please use the following code. Clone the github repository and install dependencies into a conda environment using the CompMutTB.yml file and setup Rscript. 
```
conda create -n CompMutTB
git clone https://github.com/NinaMercedes/CompMutTB.git
conda activate CompMutTB
cd CompMutTB
conda env update --file conda/CompMutTB.yml
Rscript conda/setup.R
conda deactivate
```
## How to use CompMutTB
### Pre-processing
CompMutTB is written in R language and uses custom R scripts that can be run from command line. To make CompMutTB flexible its use is split up in stages to accomodate different pre-processing pipelines. For optimal use it is recommended that CompMutTB is run using a genotype matrix. A pre-processing pipeline is available in CompMutTB an can either run from *M. tuberculosis* fastq.gz or multi-sample vcf.gz files. Please note that if you are inputting a multi-sample vcf.gz file, please ensure that your file is annotated using SnpEff (Mycobacterium_tuberculosis_h37rv) to be able to filter missense mutations, this is not a requirement if filtering by variant type is not required. An example is shown here for *katG* (missense filtering) and *oxyR'-ahpC* regions (no missense filtering). This should be done for each individual region e.g. *katG* and then *oxyR'-ahpC* as shown.
```
conda activate CompMutTB
#katG with missense filtering
Rscript data/code/pre_process.R --sample_names "data/names.txt" --vcf_gz "data/mtb_vcf/mtb.genotyped.ann.vcf.gz" --region "data/katG_regions.tsv" --out_geno "data/katG_example.geno" --drug "isoniazid" --metadata "data/test_metadata.csv" --out_geno_nolin "data/katG_nolin_example.geno" --missense TRUE --out_vcf "data/katG.vcf.gz"

#oxyR with no missense filtering
Rscript data/code/pre_process.R --sample_names "data/names.txt" --vcf_gz "data/mtb_vcf/mtb.genotyped.ann.vcf.gz" --region "data/oxyR_regions.tsv" --out_geno "data/oxyR_example.geno" --drug "isoniazid" --metadata "data/test_metadata.csv" --out_geno_nolin "data/oxyR_nolin_example.geno" --missense FALSE --out_vcf "data/oxyR.vcf.gz"

#Option to run using fastq files (not recommended)
Rscript data/code/pre_process.R --make_vcf TRUE --fastq_directory "fastq_files" --sample_names "data/names.txt" --region "data/katG_regions.tsv" --out_geno "data/katG_example.geno" --drug "isoniazid" --metadata "data/test_metadata.csv" --out_geno_nolin "data/katG_nolin_example.geno" --missense TRUE --out_vcf "data/katG.vcf.gz"
###Please note *nolin_example.geno files maybe empty due to small sample size in the example data provided, ideally should only be run with large datasets
```
The following files are required: 
- sample_names: text file containing sample names in vcf/ for fastq files
- vcf_gz: name of multi-sample vcf_gz to input
- region: tsv file regions to filter on one line only (example: Chromosome	763370	767320)
- out_geno: name of geno file to output
- drug: phenotype of interest
- metadata: csv files containing metadata for each sample, missing data should be encoded NA. Headers should include at least id,lineage,drug... (lineage not required if not filtering out lineage-specific variants --rem_lineage=FALSE)
- out_geno_nolin: name of geno file to output without lineage-specific variants
- out_vcf: name of vcf to output containing only regions data (useful for further analysis e.g. MAF calculation)

 Outputs will include a genotype matrix, filtered vcf.gz and variant information file (*_info.csv").

### Association and Mediation Analysis 
CompMutTB use both association and mediation tests to identify potential compensatory mutations. Once you have a genotype matrix for each region, you can use the following code to perform the tests (see example data for how matrices should look). Association tests are first used as a screening method to identify any mutations associated with the drug-resistant phenotype in question and their snp-snp associations. Mediation analysis follows (see paper for details on how this works). As this can take some time, it is recommended that this is only performed on significant mutation pairs output by the association analysis. A threshold can be set using the command line, if you want to run all pairs use a threshold of 1. The default is 0.05. All p-values are adjusted for using the false discovery rate.
```
#Association between rpoB and rpoC example data
Rscript data/code/assoc_test.R --geno_a "data/rpoB.geno" --geno_b "data/rpoC.geno" --chi_file "data/rpoB_rpoC_chi.csv" --metadata "data/test_metadata.csv" --drug "rifampicin" --sample_names "data/names.txt"
#Mediation between rpoB and rpoC example data with relaxed threshold of 0.6 to generate some results
Rscript data/code/med_test.R --geno_a "data/rpoB.geno" --geno_b "data/rpoC.geno" --chi_file "data/rpoB_rpoC_chi.csv" --metadata "data/test_metadata.csv" --drug "rifampicin" --sample_names "data/names.txt" --med_file "data/rpoB_rpoC_med.csv" --threshold 0.6
```
These scripts require the following files/ information:
- geno_a: a genotype matrix for potential drug-resistance mutations e.g. *rpoB* (see example data)
- geno_b: a genotype matrix for potential compensatory mutations e.g. *rpoC* (see example data)
- chi_file: name of file to output association analysis (csv)
- metadata: csv files containing metadata for each sample, missing data should be encoded NA. Headers should include at least id, drug... 
- drug: phenotype to use as outcome variable
- sample_names: text file containing sample names in metadata
- med_file: name of file to output mediation analysis (csv)
- threshold: threshold for significant association to filter out variants for mediation analysis

Outputs including results file containing results for association and mediation analyses.

### Interpretation
This step takes the final mediation analysis results and interprets how likely compensatory mutations are to be a potential compensatory mutation using pre-defined thresholds. This requires an info file (csv) containing additional context for each SNP (see example data *_info.csv), if no additional information is required, please use a file with at least two columns "snp" containing the snp ids and "info" with 'none' as a descriptor for each row. 
```
#Interpret Results
Rscript data/code/interpret_results.R   --med_file "data/rpoB_rpoC_med.csv" --out_file "data/rpoB_rpoC_example_results.csv" --chi_only "data/rpoB_rpoC_example_chi_only.csv" --info_file_r1 "data/rpoB_info.csv" --info_file_r2 "data/rpoC_info.csv"
conda deactivate
```
Inputs/ outputs include:
- med_file: name of file output by mediation analysis (csv)
- out_file: name final results file to output (csv)
- info_file_r1: variant information for potential drug-resistance mutations e.g. *rpoB* (see example data), must contain column 'snp'
- info_file_r2: variant information for potential compensatory mutations e.g. *rpoC* (see example data) (see example data), must contain column 'snp'
