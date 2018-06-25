#20180610

#sarahfong
#biopython is loaded in the sf_test env with 'conda install -c bioconda biopython'

import sys, os
import glob

#SARAH- ADD THE PATH WITH BIOPYTHON
bio_package = '/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages'
if bio_package not in sys.path:

    sys.path.append('/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages')
sys.path
    
from Bio import AlignIO
from Bio import SeqIO


# path to data hg38

datapath = '/dors/capra_lab/data_clean/ucsc/hg38/multiz100way/maf/'

outpath = '/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/hg38_human_specific_coordinates/'

#make a list of the MAF files separated by chromosome
os.chdir(datapath)
maf_list = glob.glob("chr*.maf*")
#maf_dict = {(i.split("/")[-1]).split(".")[0]: i for i in maf_list}


maf_list_short=list((i for i in maf_list if len(i)<13))
maf_dict = {(i.split(".")[0]): i for i in maf_list_short}

# make a list of villar species. Pull out any alignments where hg19 matches any of these species
#hq_species =['Mmul', 'Cjac', 'Mmus', 'Rnor', 'Ocun', 'Btau', 'Sscr', 'Cfam', 'Fcat'] #rhesus, marmoset, mouse, rat, rabbit, cow
hq_species = ['hg19', 'rheMac2', 'calJac1', 'mm9', 'rn4', 'oryCun2', 'bosTau4', 'canFam2', 'felCat3', 'susScr3']


for chr_num, MAF_FILE in maf_dict.items():
# this grafted from Abin's script. Credits to him. 
    
    # make a file for each chromosome
    outfile = "%shg38_species_counts_%s.bed" %(outpath, chr_num)

    touch = "touch %s" %outfile
    print(touch)
    os.system(touch)
    out_file = open(outfile, 'w')
    
    #unzip the maf file
    unzip_cmd = "gunzip %s"% (MAF_FILE)
    print(unzip_cmd)
    os.system(unzip_cmd)
    
    # parse the unzipped maf file
    unzipped_maf = chr_num + ".maf"
    maf = AlignIO.parse(unzipped_maf, "maf")
    count = 0 
    
    for block in maf: 
        count += 1
        store_homolog_line = []

        s_count = 0
        v_count = 0
        s_list =[]
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
            
            if species != "hg38":
                #if s_count = 0, the sequence block is human specific.
                
                s_count +=1
                s_list.append(species)
                if species in hq_species:
                    v_count +=1
            else:
                # only write the row with the h19 information. Not all the species information.
                store_homolog_line.append([blockchr, start, end, strand, species]) 
            
        store_homolog_line.append([str(s_count), str(v_count)])
        store_homolog_line = [value for v in store_homolog_line for value in v]
        store_homolog_line.append(s_list)
        
        out_file.write('\t'.join(map(str,store_homolog_line)) + '\n')
       
    out_file.close()     
    
    #rezip the maf file
    zip_cmd = "gzip %s"% (unzipped_maf)
    print(zip_cmd)
    os.system(zip_cmd)

