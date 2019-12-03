#!/bin/bash
# ------------------------------------------------------------------
# [Juan A. Rodriguez] Title run_lmi.sh
#          Description
# Runs LMI for a set of SNPs
# October 18th, 2019
# ------------------------------------------------------------------

# variables to pass:
wd=$1 # wordking directory where to send output and generate subfolders. The same where you put the summary statistics
chr_dir=$2 # Directory with all the chromosome files downloaded from 1000Genomes, filtered by CEU, and MAF >0.01.
stats=$3 # summary statistics file (column order formatted accordingly)
chr=$4 # Chromosome to run
win=$5 # Window around our SNP of interest to create a neighbourhood
mafThres=$6 # maximum difference allowed to match SNPs by MAF.
prefix=$7 # just an alias to identify your workd
exdir=$PWD # script execution directory
# Change to the working directory
cd ${wd};
# Create some folders...
mkdir -p 'sumstats';
mkdir -p 'snp_files';

#### Check your summary statistics are formatted following the order specified in the Jupyter Notebook

# Split sumstats file in chromosomes
echo -e 'Splitting in (selected) chromosomes... \n';
awk -v c=${chr} -v OFS='\t' '$2==c' ${wd}'/'${stats} | sort -gk2 > ${wd}'/sumstats/'${stats}.chr${chr}

# Run this script to generates the ecdf for the ORs
${exdir}/cdf_OddsRatio.R ${chr} ${wd} ${stats}

# Delete the LMI record file for the current chromosome
rm -f ${wd}/${stats}_localMoranI_${prefix}.chr${chr}

# Create the file to compute SNPs- Chr Pos                                                                                                  
awk -F'\t' -v OFS=$'\t' '{print $3,$1}' ${wd}'/sumstats/'${stats}.chr${chr}.cdf | sort -k1n > ${wd}'/sumstats/'${stats}.chr${chr}.snps

# START READING SNPs and calculating LMI:

echo -e 'Calculating LMI for each SNP... \n'
while read pos snp; do
    echo $snp, $pos
# Add and subtract window
    up=$(( pos + win ));
    down=$(( pos - win ));
# If the downstream position is negative, then consider 0 as a start
    if (( down < 0 )); then
        down=0;
    fi;
    echo $down , $pos , $up;
# Getting the MAF:
# Get the MAF for the viewpoint (the SNP we are looking at the moment) 
    maf_temp=$(awk -v p="${pos}" '$2==p' ${chr_dir}/ALL.chr${chr}_CEU_MAF001.frq.filt | cut -f5);
# Check MAF is less than 0.5:
# in case its higher than 0.5, then raises a flag and sets the correct value of MAF
    inv_maf=0
    read -r maf inv_maf <<< $(awk -v m=$maf_temp 'BEGIN {if(m > 0.5){print 1-m,1}else{print m,0}}')
    echo "is it inverted?:" ${inv_maf}
    echo "this is the maf:"$maf
# If the snp/position does not exist in the 1000 genomes, raise a warning and continue to the next SNP. 
    if [ -z "${maf_temp}" ]; then
        echo -e "MAF Unset!" "\n"
        continue;
    fi
# Return the MAF values for the low and up interval to match by MAF
    read -r upmaf lowmaf <<< $(awk -v m=${maf} -v mt=${mafThres} 'BEGIN{ upmaf = m + mt; lowmaf = m - mt; print upmaf; print lowmaf}')
    echo "low & up interval:"${lowmaf}, ${upmaf}
# Matching by MAF:
# For SNPs in the region from the GWAS SumStats (selecting from *.snps) do the overlap over SNPs 
# with similar MAF (selecting from *.frq.filt.), within the specified region of X Kb.
# Intersect the MAF matched selected SNPs in the region with those that have been tested in the GWAS

    awk -v OFS='\t' 'NR==FNR { a[$1]=$2; next } $2 in a { print $0,a[$2] } ' \
	<(awk -v d="${down}" -v u="${up}"  '$1 > d  && $1 < u' 'sumstats/'${stats}.chr${chr}.snps) \
	<(awk -v d="${down}" -v u="${up}" -v lf="${lowmaf}" -v uf="${upmaf}"  '$2 > d && $2 < u && $5 > lf && $5 < uf' ${chr_dir}/ALL.chr${chr}_CEU_MAF001.frq.filt) > 'snp_files/'${snp}.1KG.join;
    
# Add the values from the cdf by merging through SNP position because some of the rs's in the sumstats 
# are missing in 1KG or changed name (mostly this).
# We get the following order: CHR-POS-SNP-odds-PV-MAF
    awk -v OFS='\t' 'NR==FNR { a[$2]=$0; next } $3 in a { print $0,a[$3] } ' 'snp_files/'${snp}.1KG.join 'sumstats/'${stats}.chr${chr}.cdf | awk -v OFS='\t' '{print $2,$3,$1,$7,$5,$12}' | sort -gk2 > 'snp_files/'${snp}.1KG.sumstats;
# Generate a temporary file with the position of the neighbourhood of the SNP we are looking
    ranfile=$(mktemp /tmp/ran.XXXXXXXXXXXXX);
# Selects the corresponding SNPs for the LMI sumstats. NOW ITS SELECTING BY POSITION, field $2.
    awk -v OFS='\t' '{print $2}' 'snp_files/'${snp}.1KG.sumstats > ${ranfile};
# Calculate LD. We just need the first starting position for the SNP                                                                                        
# and the last one. 
    read first last <<< $(awk 'NR==1; END{print}' 'snp_files/'${snp}.1KG.sumstats | cut -f2);
    vcffile=$(mktemp /tmp/vcf.XXXXXXXXXXXXX);
    echo ${chr}:${first}-${last};
# using tabix, query the vcf to select the vcf with the SNPs of the neighbourhood of our SNP.
    tabix -h -f ${chr_dir}/ALL.chr${chr}_CEU_MAF001.vcf.gz ${chr}:${first}-${last} > ${vcffile};
# Add function to filter vcf and keep only the positions we are looking.
    vcffile_filt=$(mktemp /tmp/vcf.XXXXXXXXXXXXX);
    cat <(grep '#' ${vcffile}) <(awk -v OFS='\t' 'NR==FNR { a[$1]=$2; next } $2 in a { print $0 } ' ${ranfile} ${vcffile})  > ${vcffile_filt}
    plink --vcf ${vcffile_filt} --make-bed --memory 1900 --silent --out 'snp_files'/${snp};    

# For plink to calculate LD of a single variant with a set of SNPs we need to pass to plink the rs code of the variant.
# !!! WE NEED TO CHECK IF OUR SNP EXISTS WITH THAT RS_code IN THE 1KG FILE.
# OTHERWISE, THE PLINK LD COMPUTATION WILL CRASH: ->
# This piece checks the rs correct name from the 1KG vcf, sets it in variable to pass to plink.
    altsnp=$(awk -v p=$pos '{if($2==p){print $3}}' ${vcffile});
# the name of the snp we are looking
    echo 'altsnp: '${altsnp};
    plink --bfile 'snp_files/'${snp} --ld-snp ${altsnp} --ld-window-kb 10000 --ld-window-r2 0 --r2 --ld-window 1000000 --noweb --memory 1900 --silent --out 'snp_files/'${snp};
    echo -e '\n';
# Delete header and format LD output.                                                                                                                       
    awk -v OFS='\t' 'NR==FNR { a[$1]=$0; next } $2 in a { print $0,a[$2] } ' <(sed -e 1d 'snp_files/'${snp}.ld | awk -v OFS='\t' '{print $5,$6,$7}') 'snp_files/'${snp}.1KG.sumstats | cut -f1-6,9 > 'snp_files/'${snp}.1KG.sumstats.ld;
    # Calculate LMI    
    lmi=$(${exdir}/localMoranI.R ${wd}/'snp_files' ${snp});
    echo -e ${snp}"\t"${pos}"\t"${lmi} >> ${wd}/${stats}_localMoranI_${prefix}.chr${chr};
# Clean files                                                                                                                                               
    find "snp_files/" -name "${snp}.*" -delete;
    find ${vcffile} -type f -exec rm {} +;
    find ${vcffile_filt} -type f -exec rm {} +;
    find ${ranfile} -type f -exec rm {} +;
done < 'sumstats/'${stats}.chr${chr}.snps
# For testing a few SNPs...
# done < <(head -n15 'sumstats/'${stats}.chr${chr}.snps);

### END
