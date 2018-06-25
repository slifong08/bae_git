import sys, os
import glob

#SARAH- ADD THE PATH WITH BIOPYTHON
bio_package = '/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages'
if bio_package not in sys.path:

    sys.path.append('/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages')
sys.path

from Bio import AlignIO
from Bio import SeqIO

# path to data

datapath = '/dors/capra_lab/data/alignments/hg19/multiz46way/'

outpath = '/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/hg19_villar/'

#make a list of the MAF files separated by chromosome

maf_list = glob.glob("%schr*.maf"%datapath)
maf_dict = {(i.split("/")[-1]).split(".")[0]: i for i in maf_list}

# make a list of villar species. Pull out any alignments where hg19 matches any of these species
#hq_species =['Mmul', 'Cjac', 'Mmus', 'Rnor', 'Ocun', 'Btau', 'Sscr', 'Cfam', 'Fcat'] #rhesus, marmoset, mouse, rat, rabbit cow
hq_species = ['hg19', 'rheMac2', 'calJac1', 'mm9', 'rn4', 'oryCun2', 'bosTau4', 'canFam2', 'felCat3']

for chr_num, MAF_FILE in maf_dict.items():
# this grafted from Abin's script. Credits to him. 
    
    # make a file for each chromosome
    outfile = "%shg19_villar_species_v_human_%s.bed" %(outpath, chr_num)

    touch = "touch %s" %outfile
    print(touch)
    os.system(touch)
    out_file = open(outfile, 'w')
    
    maf = AlignIO.parse(MAF_FILE, "maf")
    count = 0 
    for block in maf: 
        count += 1
        store_homolog_line = []
        #print(this_paragraph)
        s_count = 0
        v_count = 0
        for row in block:             
            ### this is where the parsing happens.... what is the .id part?
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
            s_count +=1
            if species in hq_species:
                if species != "hg19":
                    v_count +=1
                else:
                    store_homolog_line.append([blockchr, start, end, strand, species])
        store_homolog_line.append([str(s_count), str(v_count)])
        store_homolog_line = [value for v in store_homolog_line for value in v]


        out_file.write('\t'.join(map(str,store_homolog_line)) + '\n')
       
    out_file.close()   