#This file parses the raw output of the focus program
#It outputs a table with the SRA ID's as headers, the Patric IDs as the first entry in each row (row names)
'''
Required inputs are
1. a directory of results files, files are named so that the ID is first and seperated by a '.'
2.the files listed in the dictionary must be tab seperated file with one header line and two or three columns. The thrid column is optional. Example format:

Strain  DRR000231_pass.fasta #header line
Firmicutes bacterium UBA5237 PATRIC|1948037.3   12.026660711310816 #line 1
Geoglobus sp. UBA232 PATRIC|1915570.3   6.136419201519487 #line 2

3. output file name 
# Note: do not run from inside the directory of results, the output will be written to the current directory
'''


import argparse
import os
import subprocess
import pandas


#open the directory and list all files
#split them to get just the SRA headers and then use this as the header line

def file_parse(path):
    dict_of_dicts={}
    for fn in os.listdir(path):
        sra_dict={}
        sra=fn.split('.')[0]
        with open(path+fn, 'r') as infile:
            next(infile)
            for line in infile:
                fields=line.strip('\n').split('\t')
                patric_id=fields[0].split('PATRIC|')[1]
                if len(fields)==2: #determine if there are one or two sets of counts in the file
                    read_count=float(fields[1])
                elif len(fields)==3: #if there is > 1 set of counts, average the two counts 
                    read_count=(float(fields[1])+float(fields[2]))/2
                elif len(fields)<2:
                    print("Improper file format")
                elif len(fields)<3:
                    print("Improper file format")
                sra_dict[patric_id]=read_count
        dict_of_dicts[sra]=sra_dict
    return(dict_of_dicts)

def file_writer(intermediate_file,data):
    df = pandas.DataFrame(data)
    df=df.fillna(0)
    df.to_csv(intermediate_file,sep='\t')
def fix_header(intermediate_file,output_file):
	i=0
	print("Writing output to " + output_file)
	with open(output_file, 'w') as outf:
		with open(intermediate_file, 'r') as inf:
			for line in inf: 
				if i==0:
					header='genome'+line
					outf.write(header)
					i+=1
				
				else:
					outf.write(line)

#write to file
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse a directory of focus files and generate a single tsv file")
    parser.add_argument('-d', help='directory of focus output files, and nothing else. Files should start with the genome ID', required=True)
    parser.add_argument('-o', help='output filename', required=True)
    args = parser.parse_args()
    intm_file='focus_file_intermediate.tsv'
    read_counts= file_parse(args.d)
    #print(sra_ids,patric_ids, read_counts)
    file_writer(intm_file,read_counts)
    fix_header(intm_file,args.o)
    subprocess.call(["rm", intm_file])
