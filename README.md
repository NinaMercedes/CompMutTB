# CompMutTB - A framework to identify compensatory mutations from *Mycobacterium tuberculosis* whole genome sequences
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
CompMutTB is written in R language and uses custom R scripts that can be run from command line. To make CompMutTB flexible its use is split up in stages to accomodate different pre-processing pipelines. For optimal use it is recommended that CompMutTB is run using a genotype matrix. A pre-processing pipeline is available in CompMutTB an can either run from *M. tuberculosis* fastq.gz or multi-sample vcf.gz files. Please note that if you are inputting a multi-sample vcf.gz file, please ensure that your file is annotated using SnpEff (Mycobacterium_tuberculosis_h37rv) to be able to filter missense mutations, this is not a requirement if filtering by variant type is not required. An example is shown here for katG (missense filtering) and OxyR'-ahpC regions (no missense filtering).
```
conda activate CompMutTB
#katG with missense filtering
Rscript data/code/pre_process.R --sample_names "data/names.txt" --vcf_gz "data/mtb_vcf/mtb.genotyped.ann.vcf.gz" --region "data/katG_regions.tsv" --out_geno "data/katG_example.geno" --drug "isoniazid" --metadata "data/test_metadata.csv" --out_geno_nolin "data/katG_nolin_example.geno" --missense TRUE --out_vcf "data/katG.vcf.gz"

#oxyR with no missense filtering
Rscript data/code/pre_process.R --sample_names "data/names.txt" --vcf_gz "data/mtb_vcf/mtb.genotyped.ann.vcf.gz" --region "data/oxyR_regions.tsv" --out_geno "data/oxyR_example.geno" --drug "isoniazid" --metadata "data/test_metadata.csv" --out_geno_nolin "data/oxyR_nolin_example.geno" --missense FALSE --out_vcf "data/oxyR.vcf.gz"

#Option to run using fastq files (not recommended)
Rscript data/code/pre_process.R --make_vcf TRUE --fastq_directory "fastq_files" --sample_names "data/names.txt" --vcf_gz "data/mtb_vcf/mtb.genotyped.ann.vcf.gz" --region "data/katG_regions.tsv" --out_geno "data/katG_example.geno" --drug "isoniazid" --metadata "data/test_metadata.csv" --out_geno_nolin "data/katG_nolin_example.geno" --missense TRUE --out_vcf "data/katG.vcf.gz"
```

