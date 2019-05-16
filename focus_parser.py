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
#add each patric ID to a set so that there is a full set of IDs, use this as the column names

#each count is stored as an entry an entry in a dictionary, with IDs as the key


#create a table with the SRA ID's as headers, the Patric IDs as the first entry in each row (row names)
#fill with zeros

#iterate through each dictionary entry, and fill in the counts in the table created above

#write to file
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse a directory of focus files and generate a single tsv file")
    parser.add_argument('-d', help='directory of focus output files, and nothing else. Files should start with the genome ID', required=True)
    parser.add_argument('-o', help='output filename', required=True)
#    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

#    test_data_path='C:\\Users\JTB\Documents\edwards_lab\projects\mgenomes_clusters\\test_data\\'
#    outf='C:\\Users\JTB\Documents\edwards_lab\projects\mgenomes_clusters\\test_data\\focus_test_data_output.tsv'
    intm_file='focus_file_intermediate.tsv'
    read_counts= file_parse(args.d)
    #print(sra_ids,patric_ids, read_counts)
    file_writer(intm_file,read_counts)
    fix_header(intm_file,args.o)
    subprocess.call(["rm", intm_file])
'''
to run from the command line:
1. provide the name of directory contain only results files
2. the results files must be named so that the SRA id is first and seperated by a '.'
3.the results files must be tab seperated file with on header line and two or three columns. The thrid column is optional. Example format:
Strain  DRR000231_pass.fasta 
Firmicutes bacterium UBA5237 PATRIC|1948037.3   12.026660711310816
Geoglobus sp. UBA232 PATRIC|1915570.3   6.136419201519487
Klebsiella pneumoniae PATRIC|573.16228  0.4643423056034224
Klebsiella pneumoniae PATRIC|573.16232  0.494672270710047
Lachnospiraceae bacterium UBA2942 PATRIC|1952019.3      5.301575914410268
Lactobacillus backii PATRIC|375175.76   0.14775219399495737
Lactobacillus delbrueckii PATRIC|767455.7       0.2672115673522007
Lactobacillus salivarius PATRIC|1624.129        0.002988559684395548
Marine Group II euryarchaeote PATRIC|2163009.67 3.19890019358266
.
.
.
4. provde the name of the output file 
#5. Note: do not run from inside the directory of results, the output will be written to the current directory
'''
