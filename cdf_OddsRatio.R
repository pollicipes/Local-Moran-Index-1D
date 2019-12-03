#!/usr/bin/Rscript

####!/apps/R/2.14.2/bin/Rscript

args<-commandArgs(TRUE)
# Chromosome number
ch <- as.character(args[1])
# The working directory
wd <- as.character(args[2])
# Summary stats name (sumstats should be placed in the working directory folder)
ss <- as.character(args[3])

# Read ss
a=read.table(paste(wd,"/sumstats/",ss,".chr",ch,sep=''))
# Get only complete cases
a.f=a[complete.cases(a),]
# Turn the OR to >1
x=ifelse(a.f$V4 < 1,1/a.f$V4,a.f$V4)
# Add vector of OR to table
a.f$sor=x
# Sort dataframe on new OR column
a.fs=a.f[order(a.f$sor),]
# Create the qvector of xth-percentiles
qvec=seq(1-0.0005,length(x)-0.0005)/length(x)
# Create the inverse cdf for the OR data
cdf=qnorm(qvec)
###hist(cdf)
# Add the data to dataframe
a.fs$cdf_or=cdf
###dim(a.fs)
# Save the dataframe to the wd/sumstats folder
write.table(a.fs, 
            file=paste(wd,"/sumstats/",ss,".chr",ch,".cdf",sep=''),
            quote=FALSE, 
            sep='\t',
            col.names = FALSE,
            row.names= FALSE)
