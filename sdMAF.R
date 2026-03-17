#!/usr/bin/env Rscript

# Copyright 2022 Zeya Chen <zeya.chen@sickkids.ca>, Zhong Wang, Delnaz Roshandel, Lei Sun, Andrew D. Paterson 
# Citation: Desmond Zeya Chen, Delnaz Roshandel, Zhong Wang, Lei Sun, Andrew D Paterson, Comprehensive whole-genome analyses of the UK Biobank reveal significant sex differences in both genotype missingness and allele frequency on the X chromosome, Human Molecular Genetics, 2023;, ddad201, https://doi.org/10.1093/hmg/ddad201
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.  
suppressPackageStartupMessages(library("argparse"))

.VERSION <- "0.0.3"

# create parser object
parser <- ArgumentParser(
  prog="sdMAF",
  description=paste0("sdMAF ",.VERSION," is a R based commend-line tool used to compute sex differences in allele frequencies. sdMAF is free and comes with ABSOLUTELY NO WARRANTY. Details of the method can be found https://journals.plos.org/plosgenetics/article/authors?id=10.1371/journal.pgen.1010231"),
  epilog="Citation: Chen et al., Comprehensive whole-genome analyses of the UK Biobank reveal significant sex differences in both genotype missingness and allele frequency on the X chromosome, Human Molecular Genetics, 2023;, ddad201, https://doi.org/10.1093/hmg/ddad201. Report bugs to zeya [dot] chen [at] sickkids [dot] ca"
  )



# specify our desired options 
# required arguments
requiredNamed = parser$add_argument_group('Required Arguments')
requiredNamed$add_argument("-f","--female", type="character", 
                           help = "Female genotype count file produced by PLINK.", 
                           metavar = "<filename>")
requiredNamed$add_argument("-m","--male", type="character", 
                           help = "Male genotype count file produced by PLINK.", 
                           metavar = "<filename>")

# optional arguments
optionalNamed = parser$add_argument_group('Optional Arguments')
optionalNamed$add_argument("--version", action="store_true",
                           help="Print the version of this tool and exit.")
optionalNamed$add_argument("--bim", type="character",
                           help = "PLINK format bim file address used to extract base pair position. Optional if ID in .gcount file are all chr:bp:a1:a2.",
                           metavar = "<filename>")
optionalNamed$add_argument("-o","--out", type="character", default="autosomal", 
                           help = "Output file name and address. Default autosomal.sdMAF in current directory. Output will look like YOURINPUT.sdMAF", metavar = "<filename>")
optionalNamed$add_argument("-l","--log", type="character",
                           help = "Log file name and address. Default 'YOURINPUTin--out'_sdMAF.log.", metavar = "<filename>")
optionalNamed$add_argument("--multi-allelic", action="store_true",default=FALSE,
                           help = "Include to keep multi-allelic SNPs in the results or not.")
optionalNamed$add_argument("--mac", type="integer", default=5,
                           help = "Sex combined minimum allele count filter. Variant with minor allele count less than input will be filtered out. Default 5.",
                           metavar = "<minimum count>")
optionalNamed$add_argument("--sex-specific", action="store_true", default=FALSE,
                           help = "Include to use sex specific minimum allele count filter for both males and females.")
  
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

#assemble log file address
if (is.null(args$log)){
  logs.nm <- paste0(args$out,"_sdMAF.log")
} else {
  logs.nm <- paste0(args$log,"_sdMAF.log")
}

#Starts logging
sink(logs.nm,split = T)

cat(paste0("########## sdMAF ",.VERSION," ########## \nAn R based commend-line tool used to compute sex differences in allele frequencies.\nsdMAF is free and comes with ABSOLUTELY NO WARRANTY.\nCitation: Desmond Zeya Chen, Delnaz Roshandel, Zhong Wang, Lei Sun, Andrew D Paterson, Comprehensive whole-genome analyses of the UK Biobank reveal significant sex differences in both genotype missingness and allele frequency on the X chromosome, Human Molecular Genetics, 2023;, ddad201, https://doi.org/10.1093/hmg/ddad201 \nReport bugs to zeya [dot] chen [at] sickkids [dot] ca.\n"))

# print version and exit early if  --version was passed
if (isTRUE(args[["version"]])){
    cat(paste0("You are using version ",.VERSION," of sdMAF tool.", "\n"))
    sink()
    quit(save = "no", status = 0)
}

cat(paste0("Analysis started at ", Sys.time(),"\n##############################", "\nChecking if inputs are valid.","\n"))
# print Error and exit early if no female genotype count found.
if (!file.exists(args[["female"]])){
  cat(paste0("Error: no female genotype count file found at ",args[["female"]],".","\n"))
  sink()
  quit(save = "no", status = 0)
}

# print Error and exit early if no male genotype count found.
if (!file.exists(args[["male"]])){
  cat(paste0("Error: no male genotype count file found at ",args[["male"]],".","\n"))
  sink()
  quit(save = "no", status = 0)
}

# loading gcount files and logging some info
cat(paste0("Loading female gcount file found from ",args[["female"]],".","\n"))
fe <- read.table(args$female)
names(fe) <- c("CHROM","ID","REF","ALT","HOM_REF_CT","HET_REF_ALT_CTS","TWO_ALT_GENO_CTS","HAP_REF_CT","HAP_ALT_CTS","MISSING_CT")
cat(paste0(sum(fe[1,5:10])," samples detected from female gcount file.","\nNumber of SNPs by chromosome from female gcount file:"))
table(fe$CHROM)
cat(paste0("Loading male genotype count file found from ",args[["male"]],".","\n"))
ma <- read.table(args$male)
names(ma) <- c("CHROM","ID","REF","ALT","HOM_REF_CT","HET_REF_ALT_CTS","TWO_ALT_GENO_CTS","HAP_REF_CT","HAP_ALT_CTS","MISSING_CT")
cat(paste0(sum(ma[1,5:10])," samples detected from male gcount file.","\nNumber of SNPs per chromosome from male gcount file:"))
table(ma$CHROM)

# print Error and exit early if female male file dimension not matching or col number is not 10.
if (nrow(fe) != nrow(ma) | ncol(fe) != 10 | ncol(fe) != 10){
  cat(paste0("Error: genotype count files dimensions are incorrect or row numbers not equal between sex.", "\n"))
  sink()
  quit(save = "no", status = 0)
}


# print Error and exit early if male genotype count file is passed to --female.
if (nrow(fe) != sum(fe$HAP_REF_CT + fe$HAP_ALT_CT == 0)){
  cat(paste0("Error: male genotype files was passed to --female.", "\n"))
  sink()
  quit(save = "no", status = 0)
}

# print Error and exit early if female male SNPs are not matching.
if (nrow(fe) != sum(fe$ID == ma$ID)){
  cat(paste0("Error: genotype count files snps are not matching.", "\n"))
  sink()
  quit(save = "no", status = 0)
}

# check the region of SNPs if they are autosomal/PAR or ChrX NPR based on input. region = 1 for autosomal/PAR and 2 for NPR. 
NPRs <- (ma$HAP_REF_CT + ma$HAP_ALT_CT != 0)&(ma$HOM_REF_CT + ma$HET_REF_ALT_CTS + ma$TWO_ALT_GENO_CTS == 0)
PARs <- (ma$HAP_REF_CT + ma$HAP_ALT_CT == 0)&(ma$HOM_REF_CT + ma$HET_REF_ALT_CTS + ma$TWO_ALT_GENO_CTS != 0)
cat(paste0(sum(NPRs)," chrX NPR SNPs detected.", "\n"))
cat(paste0(sum(PARs)," autosomal/PAR SNPs detected.", "\n"))
cat(paste0((nrow(ma)-sum(PARs)-sum(NPRs))," SNPs with zero count across all genotypes detected. Don't worry it will be filtered out later if mac > 0.", "\n"))

wald.1df.hwd.auto <- function(x)
  # 'Wald' type, 1 d.f. assuming HWD, Autosomal 
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s1+s2; r = r0+r1+r2
  pM = (0.5*s1+s2)/s; pF = (0.5*r1+r2)/r
  pAA.M = s2/s; pAA.F = r2/r
  delta.M = pAA.M-pM^2; delta.F = pAA.F-pF^2
  stat = (pM-pF)^2/(1/(2*s)*(pM*(1-pM)+delta.M)+1/(2*r)*(pF*(1-pF)+delta.F))
  -pchisq(as.numeric(stat),df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}

wald.1df.hwd.xchr <- function(x)
  # 'Wald' type, 1 d.f. assuming HWD, Xchr 
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s2; r = r0+r1+r2
  pM = s2/s; pF = (0.5*r1+r2)/r 
  pAA.F = r2/r
  delta.F = pAA.F-pF^2
  stat = (pM-pF)^2/(1/s*(pM*(1-pM))+1/(2*r)*(pF*(1-pF)+delta.F))
  -pchisq(as.numeric(stat),df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}


loop_func <- function(df){
  # df dataframe to be fed into wald.1df.hwd function.
  # reg regions 2 for autosomal/PAR and 1 for NPR.
  # pre-calculating number of snps and initialize three lists late will be used for messages.
  reg <- is.na(df[,9])
  nr <- nrow(df) 
  frac <- seq(0,1,0.05)[-1]
  prog <- paste0(as.character(frac*100),"%")
  snpc <- ceiling(frac*nr)
  fl <- c("wald.1df.hwd.xchr","wald.1df.hwd.auto") 
  LOG10P <- c()
  j = 1
  for (i in 1:nr) {
    f <- get(fl[(2-reg[i])]) #assign which function to be used based on input
    LOG10P <- c(LOG10P,f(df[i,5:10])) #appending P value
    if (i == snpc[j]) {
      if (i == nr) {cat(paste0("Finito !", "\n"))} else {cat(paste0("Now calculated ", prog[j]," (",i,"/",nr,").","\n"))} #print % of calculation done
      j = j + 1
    }
  }
  cbind(df,LOG10P)
}

# merge fe and ma to one data frame
chrom <- fe[,c("CHROM", "ID", "REF", "ALT", "HOM_REF_CT", "HET_REF_ALT_CTS", 
"TWO_ALT_GENO_CTS")]
chrom$M_A1A1.A1 <- ifelse(NPRs,ma$HAP_REF_CT,ma$HOM_REF_CT)
chrom$M_A1A2 <- ifelse(NPRs,NA,ma$HET_REF_ALT_CTS)
chrom$M_A2A2.A2 <- ifelse(NPRs,ma$HAP_ALT_CTS,ma$TWO_ALT_GENO_CTS)
names(chrom)[3:10] <- c("A1","A2","F_A1A1","F_A1A2","F_A2A2","M_A1A1.A1","M_A1A2","M_A2A2.A2")

if (is.null(args$bim)) {
  cat(paste0("No bim file provided. OK unless ID column from genotype file not all in chr:bp:A1:A2 form.","\n"))
} else {
  ch <- read.table(args$bim,header = F)
  if (nrow(ch) != nrow(chrom)) {
    cat(paste0("Error: genotype count file and bim file not having unequal number of rows.", "\n"))
    sink()
    quit(save = "no", status = 0)
  }
} 

if (!is.null(args$bim)) {
  if (ncol(ch) == 1) {
    chrom$BP <- ch$V1
  } else {chrom$BP <- ch$V4} # add BP to results from bim file
} else {chrom$BP <- sapply(strsplit(chrom$ID,":"), `[`, 2)} # get BP from ID
cat(paste0("############################## \nInput checkers all passed, now applying filters.","\n"))

#filter for only biallelic variants
if (isTRUE(args[["multi_allelic"]])){
} else {
  bia <- nchar(chrom$A2)==1&nchar(chrom$A1)==1&(!(duplicated(chrom$BP)|duplicated(chrom$BP,fromLast = T)))
  cat(paste0("Keeping ", sum(bia)," biallelic SNPs out of ", nrow(chrom)," total SNPs from Input.","\n"))
  chrom <- chrom[bia,]
}

region <- is.na(chrom$M_A1A2) #T for NPR and F for Autosomal/PAR

# getting a list of whether each SNP passes mac keeping SNPs that are 2AA + Aa and Aa + 2aa is > MAC in both sex
# for sex combined it will be 2AAf + Aaf + region*AAm + Aam and Aaf + 2aaf + Aam + region*aam is >= MAC, region is 1 for NPR and 2 for PAR.
if (isTRUE(args[["sex_specific"]])){
  macf <- ((2-region)*chrom$M_A1A1.A1+ifelse(region,0,chrom$M_A1A2) >= args$mac) & (ifelse(region,0,chrom$M_A1A2)+(2-region)*chrom$M_A2A2.A2 >= args$mac) & (2*chrom$F_A1A1+chrom$F_A1A2 >= args$mac) & (chrom$F_A1A2+2*chrom$F_A2A2 >= args$mac)
} else {macf <- ((2-region)*chrom$M_A1A1.A1+ifelse(region,0,chrom$M_A1A2)+2*chrom$F_A1A1+chrom$F_A1A2 >= args$mac) & (ifelse(region,0,chrom$M_A1A2)+(2-region)*chrom$M_A2A2.A2+chrom$F_A1A2+2*chrom$F_A2A2 >= args$mac)
}
cat(paste0("Keeping ", sum(macf)," SNPs out of ", nrow(chrom)," SNPs based on a ",ifelse(isTRUE(args[["sex_specific"]]),"sex specific","sex combined")," minor allele count filter of ",args$mac,".\n"))
chrom <- chrom[macf,]

cat(paste0("############################## \nAll filters applied, now computing sdMAF!","\n"))
# computing p value for sdMAF 
chromwithP <- loop_func(chrom)

chromwithP$Fmissing <- fe$MISSING_CT[match(chromwithP$ID,fe$ID)]
chromwithP$Mmissing <- ma$MISSING_CT[match(chromwithP$ID,ma$ID)]
chromwithP <- chromwithP[,c(1:4,11,13:14,5:10,12)] #rearrange

# compute allele frequency
region <- is.na(chromwithP$M_A1A2)
chromwithP$Ffreq <- (0.5*chromwithP$F_A1A2+chromwithP$F_A2A2)/(chromwithP$F_A1A1+chromwithP$F_A1A2+chromwithP$F_A2A2)
chromwithP$Mfreq <- (0.5*ifelse(region,0,chrom$M_A1A2)+chromwithP$M_A2A2)/(chromwithP$M_A1A1+ifelse(region,0,chrom$M_A1A2)+chromwithP$M_A2A2)
chromwithP$DIFmaf <- ifelse((chromwithP$F_A2A2+chromwithP$M_A2A2)>(chromwithP$F_A1A1+chromwithP$M_A1A1),chromwithP$Mfreq-chromwithP$Ffreq,chromwithP$Ffreq-chromwithP$Mfreq)

# some info in result
cat(paste0(nrow(chromwithP)," total SNPs in results.", "\n"))
cat(paste0(sum(region)," were chrX NPR SNPs.", "\n"))
cat(paste0(sum(!region)," were autosomal/PAR SNPs.", "\n"))
cat(paste0("Number of SNPs by chromosome table:","\n"))
table(chromwithP$CHROM)

# assemble the output file address
f.nm <- paste0(args$out,".sdMAF")

invisible(file.create(f.nm))
f <- file(f.nm, open="w") 

cat(paste0("Writing results to ",f.nm," and logs to ",logs.nm,"\n"))

write.table(chromwithP, file = f, sep = "\t", quote = F, append=FALSE, row.names = FALSE, col.names=TRUE)

cat(paste0("Analysis ended at ", Sys.time()))
sink()
close(f)
