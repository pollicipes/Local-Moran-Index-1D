# TO DO THE QQPLOT
library("qqman")

# READ THE LMI_pv SNPs, AND KEEP ONLY THOSE THAT WERE SELECTED THROUGH THE LMI THRESHOLD
a=read.table("/home/jrodriguez/scratch/credible_sets/LMI_pv.snps", header=FALSE)
colnames(a)=c("SNP","CHR","BP","OR","P","LMI")
# USE a2 AS THIS IF WE WANT THE FULL DATASET:
a2=a
# THIS RETURNS A TOTAL OF 511 SNPS, BUT 1 OF THEM HAS NO LMI, SO WE HAVE THE 510 FROM THE PAPER
a2=a[a$LMI > 5.1071,]
lmi_snps=a2[complete.cases(a2),]
#dim(lmi_snps)

# THOSE WERE ADDED BECAUSE THE PVALUE WAS SIGNIFICANT IN ANY COHORT AND WERE KEPT BECAUSE
# SURPASSED THE 1E-4 THRESHOLD HAVE A PV THAT MIGHT NOT COME FROM THE SAME COHORT THAN THE LMI.
# WE WILL ONLY CALCULATE THE CREDIBLE SETS USING THE SAME SUMMARY STATISTICS THAT WERE USED TO
# GET THE LMI CALCULATION, WHICH CORRESPOND TO THIS FILE, ALSO AVAILABLE IN THE /data/ FOLDER
# FROM THE GITHUB REPO.
sumstats=read.table("/home/jrodriguez/scratch/credible_sets/sumstats.raw.pangen_isblac_epicuro", header=TRUE)

# GET THE LIST OF SNPS TO CALCULATE CREDIBLE SETS IN A REGION OF +/-500 KB
win=500000
# CREDIBLE SET R2 THRESHOLD
R2_threshold=0.1
# CREATE A DIRECTORY
if(!dir.exists("/home/jrodriguez/scratch/credible_sets/list_vars")){
  dir.create("/home/jrodriguez/scratch/credible_sets/list_vars", recursive = TRUE)
}

# DEFINE THE FUNCTION TO GET THE CREDIBLE SETS.
# COPIED FROM https://github.com/hailianghuang/FM-summary
getCredibleSNP <- function(snp, logProb, threshold=0.99){
  
  prob <- exp(logProb) 
  prob_normed <- prob/sum(prob)
  
  prob_cumsum <- cumsum(sort(prob_normed, decreasing=TRUE))
  nSNP <- as.numeric(which.max(prob_cumsum>threshold ))
  
  if(nSNP<=0){ nSNP=1 }
  
  credible_set <- snp[order(prob_normed, decreasing=TRUE)[1:nSNP]]
  select <- rep(FALSE, length(snp))
  select[order(prob_normed, decreasing=TRUE)[1:nSNP]] <- T 
  ret <- list(nSNP = nSNP, 
              prob_normed = prob_normed, 
              prob_cumsum = prob_cumsum[rank(-prob_normed, ties.method="random")], 
              credible_set=credible_set, 
              select=select 
  )
}
std <- function(x) sd(x)/sqrt(length(x))

# SET AN OBJECT TO STORE THOSE SNPS THAT HAVE BEEN ALREADY PICKED IN A CREDIBLE SET.
picked=c()
# HOW MANY UNIQUE SNPS WE WILL CALCULATE
total=0
# IN HOW MANY CASES ANY SNP FROM THE CREDIBLE SET IS ALSO THE TOP PVALUE SNP
count=0
# EMPTY DATAFRAME TO STORE INFO FOR THE SNPS
aa=data.frame()
# START READING THE SNPS
for (i in seq(1,dim(lmi_snps)[1])){
  # GET THE CORRESPONDING DATA
  s=lmi_snps[i,]
  # IF THE SNP HAS ALREADY BEEN PICKED IN A CREDIBLE SET PREVIOUSLY THEN SKIP IT.
  if (as.character(s$SNP) %in% picked){
    cat(as.character(s$SNP),"ALREADY PICKED!\n")
    next;
  }
  total=total+1
  # GET THE LIST OF POSITIONS IN THE 500kb +/- WINDOW TO CALCULATE LD
  up=s$BP+win
  down=s$BP-win
  # GET THE LIST OF SNPs FROM THE SUMSTATS
  list_vars=sumstats[sumstats$CHR == s$CHR & sumstats$BP >= down & sumstats$BP <= up,]
  # ORDER THEM
  list_vars=list_vars[order(list_vars$BP),]
  # GET THE ONES WITH NO MISSING VALUES
  list_vars=list_vars[complete.cases(list_vars),]
  # MAKE A BED WITH THE POSITIONS TO CALCULATE LD
  vars_bed=data.frame(list_vars$CHR,list_vars$BP-1,list_vars$BP,"0")
  vars_bed_path=paste0("/home/jrodriguez/scratch/credible_sets/list_vars/list_vars_",s$CHR,"_",s$BP,"_",s$SNP,".bed")
  write.table(vars_bed,file=vars_bed_path,quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
  # THE FOLDER WITH THE PLINK BED FILES FOR THE SNPS, AS USED IN THE PAPER
  chr_dir="/home/jrodriguez/scratch/chr_files_1KG_CEU_MAF001/"
  # PATH TO STORE RESULTS
  plink_bed_out=paste0("/home/jrodriguez/scratch/credible_sets/list_vars/plink_",s$CHR,"_",s$BP,"_",s$SNP)
  # PLINK CALL TO CREATE THE BED FILE WITH ONLY THE VARIANTS IN THE REGION
  try(system(paste0("/home/jrodriguez/plink2 --bfile ",chr_dir,"ALL.chr",s$CHR,"_CEU_MAF001 --extract range ",vars_bed_path," --make-bed --silent --out ", plink_bed_out), intern = TRUE, ignore.stderr = FALSE))
  # GET THE NAME (RS CODE) OF THE VARIANT IN CASE IT IS NOT IN THE ARRAY AS SUCH
  greped=try(system(paste0("grep ", s$BP," ",plink_bed_out,".bim"), intern = TRUE, ignore.stderr = TRUE))
  rs=unlist(strsplit(greped,"\t"))[2]
  # CALCULATE LD WITH THE SELECTED SNP
  try(system(paste0("/home/jrodriguez/plink_linux_x86_64/plink  --bfile ", plink_bed_out," --r2 --ld-snp ",rs," --ld-window-kb ",win/1000, " --silent --ld-window 99999 --ld-window-r2 0 --out ", plink_bed_out), intern = TRUE, ignore.stderr = FALSE))
  ld=read.table(paste0(plink_bed_out,".ld"),header=TRUE)
  dat=merge(list_vars,ld[,c(5,6,7)], by.x=c("BP"),by.y=c("BP_B"))
  # CALCULATE CREDIBLE SETS (Adapted from the original script for credible sets)
  select=!is.na(dat$R2) & dat$R2 > R2_threshold
  #tag_snp_P=s$P
  #tag_snp=s$SNP
  ret=getCredibleSNP(as.character(dat$SNP[select]), qchisq(dat$P[select], 1, low=F)/2)
  inCredible=rep(NA, length(dat$P))
  inCredible[match( dat$SNP[select], dat$SNP )]=0
  inCredible[match( dat$SNP[select][ret$select], dat$SNP )]=1
  prob_norm=rep(NA, length(dat$P))
  prob_norm[match(dat$SNP[select], dat$SNP )]=ret$prob_normed
  prob_cumsum=rep(NA, length(dat$P))
  prob_cumsum[match(dat$SNP[select], dat$SNP )]=ret$prob_cumsum
  result=cbind(dat, inCredible=inCredible, probNorm=prob_norm, cumSum=prob_cumsum)
  result=result[,c(3,1,2,6,4,5,7,8,9,10)]
  colnames(result)=c("CHR","BP","array_SNP_ID","SNP","OR","P","R2","inCredible","probNorm","cumSum")
  # CHECK THAT THE SEED VARIANTS WE PICK IN A CREDIBLE SET ARE NOT USED AS SEED FOR THE NEXT CREDIBLE SETS
  # WE WILL GENERATE A VECTOR CONTAINING THE SELECTED SNPS TO INCLUDE IN THE CREDIBLE SET.
  pseeds=as.character(result[which(result$inCredible == 1),]$array_SNP_ID)
  picked=c(picked,pseeds)
  # CHECK HOW MANY OF THE SNPS ARE ALSO THE TOP PV IN THE REGION
  toppv=result[which.min(result$P),]$inCredible
  if (!is.na(toppv)){
    count=count+1
  }
  # WE WILL GENERATE NOW THE DATAFRAME SUMMARIZING THE CREDIBLE SETS.
  # IT WILL INCLUDE:
  # BEFORE, DO SELECT ONLY THE SNPS CONTAINED IN CREDIBLE SET
  cset=result[which(result$inCredible == 1),]
  # TOP SNP FROM CREDIBLE SET
  ts=cset[which.min(cset$P),][,c(1:6)]
  # IS IT THE TOP IN THE REGION?
  ts$is_top=toppv
  # VARIANTS IN CREDIBLE SET
  vars_in_cset=dim(cset)[1]
  ts$vars=vars_in_cset
  # MEDIAN OR
  ts$med_OR=median(ifelse(cset$OR < 1,1/cset$OR, cset$OR))
  # STANDARD ERROR OF THE OR.
  ts$se=std(ifelse(cset$OR < 1,1/cset$OR, cset$OR))
  aa=rbind(aa,ts)
  # SAVE THE DATA FRAME WITH THE NAME REFERRING TO THE LOWEST PV SNP IN THE SET
  write.table(result, file=paste0("/home/jrodriguez/scratch/credible_sets/results/",s$CHR,"_",s$BP,"_",s$SNP,"_credible.txt"), sep="\t", quote=F, col.names=T, row.names=F)
}

### QQPLOTS FOR THE 510 SNPS CATCHED THROUGH LMI. ###

# WE PREVIOUSLY CAT THEM ALL TOGETHER, IN BASH:
#$> cat *credible.txt | grep -v 'CHR' > all_199_credible.txt
# READ TABLE
credible_regions=read.table("/home/jrodriguez/scratch/credible_sets/results/only_510_lmi/all_199_credible.txt")
colnames(credible_regions)=c("CHR","BP","array_SNP_ID","SNP",	"OR",	"P",	"R2",	"inCredible",	"probNorm", "cumSum")
head(credible_regions)
# QQPLOT FOR THE FULL REGION
pdf("/home/jrodriguez/scratch/credible_sets/qqplot_full_1Mb_region_510_SNPs.pdf")
qq(credible_regions$P,pch=20,cex=0.8,col="blue",main="Full region: LMI SNP +/- 500kb",xlim=c(0,7),ylim=c(0,7))
dev.off()
# QQPLOT ONLY FOR THE SNPS IN CREDIBLE SETS
credible_sets=credible_regions[which(credible_regions$inCredible == 1),]
dim(credible_sets)
pdf("/home/jrodriguez/scratch/credible_sets/qqplot_only_credible_510_SNPs.pdf")
qq(credible_sets$P,pch=20,cex=0.8,col="red",main="Credible sets from the LMI SNP",xlim=c(0,7),ylim=c(0,7))
dev.off()

### HOW MANY TIMES THE MOST SIGNIFICANT SNP IN THE REGION IS IN THE CREDIBLE SET

dim(aa[which(aa$is_top == 1),])
183/281
dim(aa)

###
# FOR THE FUMA ANNOTATION:
# DIVIDE INTO THE SIGNIFICANT AND NON SIGNIFICANT LEAD SNPS
cla="all_624_pv_lmi"
nosig=aa[aa$P >= 1e-4,]
write.table(nosig,file=paste0("/home/jrodriguez/scratch/credible_sets/results/",cla,"/nosig_1e4.bed"),sep="\t",quote=FALSE,row.names = FALSE)
# FORMAT FOR LEAD SNP IN FUMA
lead_nosig=nosig[,c(4,1,2)]
write.table(lead_nosig,file=paste0("/home/jrodriguez/scratch/credible_sets/results/",cla,"/lead_nosig_1e4.bed"),sep="\t",quote=FALSE,row.names = FALSE)

sig=aa[aa$P < 1e-4,]
write.table(sig,file=paste0("/home/jrodriguez/scratch/credible_sets/results/",cla,"/sig_1e4.bed"),sep="\t",quote=FALSE,row.names = FALSE)
# FORMAT FOR LEAD SNP IN FUMA
lead_sig=sig[,c(4,1,2)]
write.table(lead_sig,file=paste0("/home/jrodriguez/scratch/credible_sets/results/",cla,"/lead_sig_1e4.bed"),sep="\t",quote=FALSE,row.names = FALSE)

# THIS WILL CREATE THE LEAD SNPS FILE TO UPLOAD TO FUMA.
# (ITS JUST THE RESULT OF MERGING THE SIGNIFICANT WITH THE NON SIGNIFICANT AT 10e-4, BUT IT IS SEPARATED IN TWO FILES.)
#$> cat lead_sig_1e4.bed lead_nosig_1e4.bed | grep -v 'SNP' > lead_all.txt
#$> cat *credible* | grep -v 'CHR' |  awk '$8 == 1' > all_only_in_credible_624.txt

### END ###

# THIS CODE BELOW GETS THE CREDIBLE SETS FOR THE SNPs THAT WERE SIGNIFICANT OR NOT.
# WE WOULD NEED TO CHANGE THE NAME OF THE *_credible.txt
# BECAUSE THE SNP IN THE NAME, WHICH IS THE LMI SELECTED ONE
# MIGHT NOT CORRESPOND TO THE ONE WITH THE MINIMUM PVALUE IN THE 
# CREDIBLE SET, WHICH IS THE ONE ANNOTATED IN aa.

# fsig=paste0(sig$CHR,"_",sig$BP)
# # LIST OF FILES
# lfsig=c()
# # ITERATE THE LIST OF FILES.
# for (i in fsig){
#  lfsig=c(lfsig,list.files("/home/jrodriguez/scratch/credible_sets/results/", pattern=i, recursive = FALSE, full.names = TRUE))
# }
# mybeds_sig=do.call("rbind", lapply(lfsig, read.delim))
# mybeds_sig=mybeds_sig[which(mybeds_sig$inCredible == 1),]
# dim(mybeds_sig)
# write.table(mybeds_sig,file="/home/jrodriguez/scratch/credible_sets/results/credible_sig_1e4.bed",sep="\t",quote=FALSE,row.names = FALSE)
# 
# # CODES CHR_POS TO GET THE CURRENT FILES
# fnosig=paste0(nosig$CHR,"_",nosig$BP)
# # LIST OF FILES
# lfnosig=c()
# # ITERATE THE LIST OF FILES.
# for (i in fnosig){
#  lfnosig=c(lfnosig,list.files("/home/jrodriguez/scratch/credible_sets/results/", pattern=i, recursive = FALSE, full.names = TRUE))
# }
# mybeds_nosig=do.call("rbind", lapply(lfnosig, read.delim))
# mybeds_nosig=mybeds_nosig[which(mybeds_nosig$inCredible == 1),]
# write.table(mybeds_nosig,file="/home/jrodriguez/scratch/credible_sets/results/credible_nosig_1e4.bed",sep="\t",quote=FALSE,row.names = FALSE)
