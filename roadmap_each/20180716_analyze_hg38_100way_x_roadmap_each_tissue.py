#180716

#!/usr/bin/python

#sarahfong

# The purpose of this script is to splite each roadmap sample enhancer dataset into chromosomes 
#and intersect them against Hg38 100way MultiZ synteny blocks.

import sys, os
import subprocess
import glob
import timeit


# load the bedtools module
#os.system("ml restore bedR_ml")



# this directory contains the .bed files of all human-specific coordinates of the genome. 
hg38_dir = "/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/hg38_human_specific_coordinates/"



# get hg19 human specific bed files
hg38_list = glob.glob("%shg38_species*.bed" % hg38_dir)
hg38_dict = {(((i.split("/")[-1])).split("_")[-1]).split(".")[0]: i for i in hg38_list} # {chr:/file/path/chr.bed}
print(hg38_dict.keys())


# get roadmap functional enhancers bed files
bedpath = '/dors/capra_lab/users/fongsl/roadmap/hg38_2018-06-20/'
os.chdir(bedpath)


bed_list = glob.glob("%sHsap_H3K27ac_plus_H3K4me3_minus_*.bed"% bedpath)

bed_dict = {(i.split("/")[-1]).split("_")[-1]: i for i in bed_list} # {chr:/file/path/chr.bed}

# make a function to split bedfiles by chromosome for intersection with synteny blocks.
def synteny_query(sample, bedfile):

    start = timeit.default_timer()
    
    # split up the bedfile by chromosome to intersect with synteny blocks.
    awk_cmd='awk -F, \'{print > $1\".bed\"}\' FS="\\t" %s' % bedfile#.split("/")[-1]
    subprocess.call(awk_cmd, shell=True)

    #concatenate the chromosome files we just split up. 
    chr_files=glob.glob("%schr*" %bedpath)
    sample = sample.split(".")[0]
    print("working on sample ", sample)
    
    for b_file in chr_files: 
        # prepare to cross the enhancer dataset with the 100way multiZ synteny block
        b_chr = (b_file.split("/")[-1]).split(".")[0]

        if b_chr in hg38_dict.keys():

            # find the synteny block BED file to compare. 
            hg38_file = hg38_dict[b_chr]

            #make an intersection file
            outfile = "%s%s_%s_x_100way.bed" %(bedpath, b_chr, sample)
    
            # BED intersect where every single base-pair overlap is accounted for.
            bed_cmd = "bedtools intersect -a %s -b %s -wao > %s" %(b_file, hg38_file, outfile)
            #print(bed_cmd)
            os.system(bed_cmd)

        # toss the alternative and unalignable regions of the human genome
        clean_up = "rm %s" % b_file
        #print(clean_up)
        os.system(clean_up)

    #concatenate all of the chromosomes back together. 
    concat_cmd = "cat %schr*_%s_x_100way.bed > %senh_x_hg38_100way/%s_x_100way.bed"% (bedpath, sample, bedpath, sample)
    print(concat_cmd)
    os.system(concat_cmd)

    #remove the chromosome files for a single bedfile
    rm_cmd = "rm %schr*.bed" % bedpath
    print(rm_cmd)
    os.system(rm_cmd)
    
    stop = timeit.default_timer()
    print("it took this many seconds to process sample ", stop-start)


for sample, bedfile in bed_dict.items():
    synteny_query(sample, bedfile)

