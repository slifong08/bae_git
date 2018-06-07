
#!/bin/python
#180601
#sarahfong
#The purpose of this script is to generate a way to query broadly active enhancer BED files 
#against multiple sequence alignments. 

#Multiple Sequence Alignment files: I will use 46-way UCSC hg19 MAFs from /dors/capra_lab/data_clean/ucsc/hg19/maf/


import os
import glob
import pandas

mafpath ='/dors/capra_lab/data/alignments/hg19/multiz46way/'

#human-specific villar reference by chr number 
#chr files generated with /dors/capra_lab/users/fongsl/broadly_active_enhancers/bae_git/roadmap_multi/separate_short_multi_intersect_by_chr.ipynb 
datapath = "/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/villar_ref_data/"

## LOAD PHAST
os.system("ml load  GCC/5.4.0-2.26"
os.system("ml load Intel/2016.3.210-GCC-5.4.0-2.26")
os.system("ml load PHAST/1.4")

# make a list of the testbed files

testbed_list = glob.glob("%schr*"%datapath)

# get rid of the x-chromosome.
testbed_list.remove("/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/villar_ref_data/chrX_hspecific_villar_bae_hq_aln.bed")

print((testbed_list))

# gather the MAF files
maf_list = glob.glob("%schr*.maf"%mafpath)
print(maf_list)


maf_dict = {(i.split("/")[-1]).split(".")[0]: i for i in maf_list}
print(maf_dict)


#chr_number = 'chr1', 'chr11', etc.
def unzipzip(chr_num, filename):
    
    parse_cmd= "maf_parse -g %s %s > %sparsed_%s.maf"%(testbed, filename , datapath, chr_num)
    os.system(parse_cmd)

for testbed in testbed_list:
    chr_num= (testbed.split('/')[-1]).split('_')[0]
    filename = maf_dict[chr_num]
    print(chr_num)
    unzipzip(chr_num, filename)

