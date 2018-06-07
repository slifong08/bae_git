
# coding: utf-8

# In[1]:


#180601

#sarahfong
#biopython is loaded in the sf_test env with 'conda install -c bioconda biopython'


import sys
import glob

#SARAH- ADD THE PATH WITH BIOPYTHON
bio_package = '/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages'
if bio_package not in sys.path:

    sys.path.append('/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages')
sys.path

from Bio import AlignIO
from Bio import SeqIO

# path to data

datapath = "/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/villar_ref_data/villar_maf/"

#make a list of the MAF files separated by chromosome

maf_list = glob.glob("%sparsed*.maf" % datapath)


for MAF_FILE in maf_list:
# this grafted from Abin's script. Credits to him. 

    chr_num = ((MAF_FILE.split('/')[-1]).split('_')[1]).split('.')[0]
    print("working on ", chr_num)
    out_file = open("%shg19_spec_%s.bed" %(datapath, chr_num), 'w')

    maf = AlignIO.parse(MAF_FILE, "maf")
    count = 0 
    for block in maf: 
        count += 1
        store_homolog_line = []
        conservation_count = len(block)
        for row in block: 
                
            ### this is where the parsing happens.... what is the .id part?
            #print(row)
            species = row.id.split('.')[0]
        
            blockchr = row.id.split('.')[1]
        
            start = row.annotations['start']
            
            end =  row.annotations['start'] + row.annotations['size'] + 1 
            
            if row.annotations['strand'] == 1:
                strand = "+" 
            elif row.annotations['strand'] == -1: 
                strand = "-"
            else: 
                raise ValueError('strand parsing did not work')     
            

            store_homolog_line.append([blockchr, start, end, strand, species, conservation_count])
            
            store_homolog_line = [value for v in store_homolog_line for value in v]
            out_file.write('\t'.join(map(str,store_homolog_line)) + '\n')

