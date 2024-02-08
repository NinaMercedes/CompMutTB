## Pre-process Data
# Load required packages 
library(optparse)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
# Set seed
set.seed(10)
option_list <- list( 
  make_option(c("-z", "--make_vcf"), type="character", default=FALSE,
    help="Make multi-vcf file from fastq: arguments should be TRUE or FALSE"),
  make_option(c("-q", "--fastq_directory"), type="character", default=".",
    help="Directory containing fastq-files"),
  make_option(c("-s", "--sample_names"), type="character", default="sample_names.txt",
    help="Name of file containing sample names"),
  make_option(c("-v", "--vcf_gz"), type="character", default="mtb_vcf/mtb.genotyped.ann.vcf.gz",
    help="Insert path to vcf.gz file"),
  make_option(c("-r", "--region"), type="character", default="region.tsv",
    help="Insert path to tsv file containing gene boundaries"),
  make_option(c("-f", "--filter_vcf"), type="character", default=TRUE,
    help="Filter vcf file: arguments should be TRUE or FALSE"),
  make_option(c("-o", "--out_vcf"), type="character", default="region.vcf.gz",
    help="Name of filtered vcf.gz to output for gene "), 
  make_option(c("-j", "--make_geno"), type="character", default=TRUE,
    help="Make genotype file: arguments should be TRUE or FALSE"), 
  make_option(c("-g", "--out_geno"), type="character", default="geno_a.geno",
    help="Name of geno file to output"), 
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",
    help="File containing metadata with DST information"),
  make_option(c("-d", "--drug"), type="character", default="rifampicin",
    help="Column name for drug DST information in metadata file"), 
  make_option(c("-l", "--rem_lineage"), type="character", default=TRUE,
    help="Remove lineage SNPs: arguments should be TRUE or FALSE"), 
  make_option(c("-n", "--out_geno_nolin"), type="character", default="geno_a_nolin.geno",
    help="Name of geno file to output with no lineage SNPs"),
  make_option(c("-x", "--missense"), type="character", default=TRUE,
    help="Filter missense SNPs: arguments should be TRUE or FALSE"))
# Parse Options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
#####ADD FUNCTIONS HERE
# Function to make multi-sample vcf and annotate
make_vcf <- function(sample_names, directory){
cat(" *Making vcf \n")
  sample_name <- read.table(sample_names)
  sample_name <- sample_name[,1]
    system(paste0("mkdir mtb_vcf"))
  for (f in sample_name){
    system(paste0("fastq2vcf.py all --read1 ",directory,"/",f,"_1.fastq.gz --read2 ",directory,"/",f,"_2.fastq.gz --ref H37Rv.fasta --prefix mtb_vcf/",f," --redo")) #get vcf files 
  }
  system(paste0("mkdir mtb_vcf"))
  system(paste0("merge_vcfs.py import --sample-file ",sample_names," --ref H37Rv.fasta --prefix mtb_vcf/mtb  --vcf-dir mtb_vcf")) #make multi-vcf
  system(paste0("merge_vcfs.py genotype --ref H37Rv.fasta --prefix mtb_vcf/mtb")) #set genotypes
  system(paste0("mv mtb_vcf/*.genotyped.vcf.gz mtb_vcf/mtb.genotyped.vcf.gz")) #rename file, remove data
  system(paste0("filter_tb_vcf.py --vcf mtb_vcf/mtb.genotyped.vcf.gz --ref H37Rv.fasta")) #filter tb and fix genotypes using allelic depth
  system(paste0("bcftools annotate --rename-chrs rename_chr.txt -Oz mtb_vcf/mtb.genotyped.filtered_no_indels.vcf.gz > mtb_vcf/mtb.genotyped.chr.vcf.gz")) #rename chromosomes to be compatible with snpeff
  system(paste0("snpEff download -v Mycobacterium_tuberculosis_h37rv")) #download snpeff db
  system(paste0("snpEff ann -v Mycobacterium_tuberculosis_h37rv mtb_vcf/mtb.genotyped.chr.vcf.gz > mtb_vcf/mtb.genotyped.ann.vcf")) #annotate using snpeff
  system(paste0("bgzip mtb_vcf/mtb.genotyped.ann.vcf")) #rezip vcf
  #system(paste0("rm -r mtb*")) #remove directories that are no longer required
}
# Function to filter vcf by gene region and get missense SNPs
filter_vcf <- function(region, vcf_gz, out_vcf, missense){
  cat(" *Filtering vcf \n")
  if(missense==TRUE){
    arg1 <- capture.output(cat("SnpSift filter \"ANN[0].EFFECT has 'missense_variant'\""))
    arg2 <- vcf_gz
    arg3 <- "> temp.vcf"
    system(paste(arg1, arg2, arg3))
    system(paste("bgzip -c temp.vcf > temp.vcf.gz"))
    system(paste("rm temp.vcf"))
    vcf_gz <- "temp.vcf.gz"
  } 
  system(paste("bcftools index -c --force", vcf_gz))
  system(paste("bcftools view -Oz -R", region, vcf_gz, " > regions.vcf.gz"))
  system(paste("bcftools index -c --force regions.vcf.gz"))
  system(paste("bcftools norm -Oz -m-any regions.vcf.gz >", out_vcf))
  if(missense==TRUE){
    system(paste("rm temp.vcf.gz"))
  }
  system(paste("rm regions.vcf.gz"))
  make_info(out_vcf)
}
# Function to make info file
make_info <- function(out_vcf){
  system(paste0("bcftools query -f '%CHROM\t%POS\t%ALT\t%ANN\n' ",out_vcf," > mtb.annotation.txt"))
  anno <- read.table("mtb.annotation.txt")
  anno <-within(anno, V4<-data.frame(do.call('rbind', strsplit(as.character(V4), ',', fixed=TRUE))))
  anno <- data.frame(anno$V1, anno$V2, anno$V3, anno$V4$X1)
  colnames(anno) <- c("chr", "pos", "alt", "anno")
  anno <- anno %>% separate(anno, sep="[|]", c("alt2", "type", "effect", "gene_name", "gene", "transcript", "transcript_id", "protein_coding","PC_N","NA_change", "AA_change", "N1", "N2", "N3","SP"))
  anno <- anno[,-4]
  anno <- anno %>% select("chr","pos","alt", "type", "effect","gene_name","gene","protein_coding", "NA_change", "AA_change" )
  anno$snp <- paste0(anno$chr,"_",anno$pos,"_",anno$alt)
  fn <- gsub(".vcf.gz", "", out_vcf)
  write.csv(anno,paste0(fn,"_info.csv"), row.names=FALSE)
}
# Function to make genotype file
make_geno <- function(out_vcf, geno, sample_names){
  cat(" *Making geno file  \n")
  system(paste("bcftools index -c --force",out_vcf))
  system(paste("bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]\n'", out_vcf, " | tr '|' '/' | sed 's/\\.\\/\\./Ns/g' | sed 's/0\\/1/0.5/g' | sed 's/[123456789]\\/[123456789]/1/g' | sed 's/0\\/0/0/g' >", geno))
  system(paste("bcftools query -l", out_vcf, " >", sample_names))
  system(paste("awk \'BEGIN{OFS=\"_\"} {print $1,$2,$3}\'", geno, "> tmp"))
  system(paste("awk \'NR==FNR{a[NR]=$0;next}{$1=a[FNR]}1\' tmp ", geno, ">tmp2 &&  cut --complement -d' ' -f 2,3 tmp2 >tmp3 && mv tmp3", geno, "&& rm tmp2 && rm tmp"))  
}
# Function to format genotype
format_geno <- function(geno_table, sample_names, out_geno){
  cat(" *Formatting genotypic data \n")
  geno_table <- read.table(geno_table)
  sample_name <- read.table(sample_names)
  snp_name <- geno_table[,1]
  geno_table <- geno_table[,-c(1)]
  is.na(geno_table)<- geno_table == "Ns" 
  geno_table <- t(geno_table)
  colnames(geno_table)<- snp_name
  geno_table <- apply(geno_table,2, function(x) as.numeric(as.character(x)))
  columns <- names(which(colSums(geno_table, na.rm=TRUE)>0))
  geno_table <- geno_table[,(which(colSums(geno_table, na.rm=TRUE)>0))]
  geno_table <- data.frame(geno_table)
  colnames(geno_table) <- columns
  rownames(geno_table) <- sample_name[,1]
  write.table(geno_table, out_geno, sep = " ", quote=FALSE, row.names=TRUE)
  return(geno_table)
}
# Functions to remove lineage SNPs
lineage_loop <- function(geno_table, metadata, drug){
  metadata <- metadata[grep("1", metadata[,drug]),]
  lineage <- matrix(ncol=1, nrow=ncol(geno_table))
  lineage <- data.frame(lineage)
  mylist <- list()
  for (i in colnames(geno_table)){
    mylist[[i]] <- geno_table[grep("1",geno_table[,i]),]
    mylist[[i]] <- rownames(mylist[[i]])
    mylist[[i]] <- metadata[metadata$id %in% mylist[[i]],]
    mylist[[i]]<- nrow(unique(mylist[[i]]['lineage']))
    if (mylist[i] < 2 ){lineage[i,1] <- names(mylist[i])
    } 
  }
 return(lineage)
}
# Function to remove lineage-specific SNPs
rem_lin <- function(metadata, drug, geno, geno_nolin_file, sample_names, out_geno){
  metadata <- read.csv(metadata, header=TRUE)
  #metadata$lineage <- gsub("lineage1.*", "lineage1", metadata$lineage)
  #metadata$lineage <- gsub("lineage2.*", "lineage2", metadata$lineage)
  #metadata$lineage <- gsub("lineage3.*", "lineage3", metadata$lineage)
  #metadata$lineage <- gsub("lineage4.*", "lineage4", metadata$lineage)
  #metadata$lineage <- gsub("lineage5.*", "lineage5", metadata$lineage)
  #metadata$lineage <- gsub("lineage6.*", "lineage6", metadata$lineage)
  #metadata$lineage <- gsub("lineage7.*", "lineage7", metadata$lineage)
  geno <- format_geno(geno, sample_names, out_geno)
  lineage <- lineage_loop(geno, metadata, drug)
  lineage <- data.frame(lineage)
  colnames(lineage) <- c("lin")
  lineage <- lineage[!is.na(lineage$lin),]
  geno_no_lin <- geno[,-which(colnames(geno) %in% lineage)]
  write.table(geno_no_lin, geno_nolin_file, sep = " ", quote=FALSE, row.names=TRUE)
}
# Run
if(opt$make_vcf==TRUE){
  make_vcf(opt$sample_names, opt$fastq_directory)
}
if(opt$filter_vcf==TRUE){
  filter_vcf(opt$region, opt$vcf_gz, opt$out_vcf, opt$missense)
}
if(opt$make_geno==TRUE){
  make_geno(opt$out_vcf, opt$out_geno, opt$sample_names )
}
if(opt$rem_lin==TRUE){
  rem_lin(opt$metadata, opt$drug, opt$out_geno, opt$out_geno_nolin, opt$sample_names, opt$out_geno)
}