#This script contains functions the following functions:

# 1. get_ids - parse a csv file that has two columns: 1) "ground truth" classifications & 2) a list of SRA runs (comma seperated)
# 2. get_focus - parse the PARTIE_HAT (FOCUS) output (.json) and obtain the FOCUS data for each ID from part 1
# 3. get_partie - parse the the PARTIE output (.tsv) to obtain entries for each ID from part 1
# 4. format_data - format all the data from part 1-3 into a table (.csv) order to be imported into R for further analysis

# it takes as input arguments:
# 1. mininum number of runs (this is the number of runs from each class to include in the output
# 2. 3 file paths - (a) to the data required by function 1, "get_ids" (b) to the data required by function 2, the output of FOCUS & (c) to the data required by function 3, the output of PARTIE 
# 3. an output file path

#it outputs a single file that contains a balanced data set of all input features from partie and focus that can be used in subsequent analysis
import argparse
import csv
import random 
# 1.
def get_ids(id_file):
 run_ids={}
 with open(id_file,'r', encoding='utf-8') as infile:
  next(infile)
  for line in infile:
    line=line.split(',')
    env = line[0].strip(' ')
    ids=line[1:len(line)]
    ids1= [x.strip('\n').strip('"') for x in ids]
    if env in run_ids.keys():
     for r in range(0,len(ids1)):
      run_ids[env].append(ids1[r])
    else:
      run_ids[env]=list()
      for r in range(0,len(ids1)):
       run_ids[env].append(ids1[r])
 for env in run_ids:
     all_ids=run_ids[env]
     run_ids[env]=random.sample(all_ids, len(all_ids)) #this function places the runs IDs from "get_ids" function in a randomized order
 return(run_ids)    

# 2.
#for each SRA run ID in the dictionary of id (from part 1):
#create get the row of info from the focus tsv and store it in a new dictionary with the ID as the key
#the arugment "minimum" includes onlu environments with that at least that many of IDs associated with it
#and it includes only that number of IDs. So if the minimum is set to 50 is will exclude environments with less than 50 entries
#and it will include only 50 entries for those environments that meet the cutoff

def get_focus(run_ids, minimum, focus_runs):
 print("Parsing focus output...")
 sra_ids = list()
 focus_nums={}
 with open(focus_runs,'r') as focus_file:
  for env in run_ids:
   if len(run_ids[env]) > minimum-1: 
    i=0
    if minimum==1:
      for sra_id in run_ids[env]:
          sra_ids.append(sra_id)
          focus_file.seek(0)
          for line in focus_file: 
           if sra_id in line:
            line=line.split('\t')
            sra_id_json= line[0]
            entries = line[1].strip('\n').strip('"')
            focus_nums[sra_id]=entries
    else:
        while i<minimum:
            sra_id=run_ids[env][i]
            sra_ids.append(sra_id)
            focus_file.seek(0)
            i+=1
            for line in focus_file: 
               if sra_id in line:
                line=line.split('\t')
                sra_id_json= line[0]
                entries = line[1].strip('\n').strip('"')
                focus_nums[sra_id]=entries

 sra_ids=set(sra_ids)       
 print(len(focus_nums.keys()), ' out of', len(sra_ids),' matching entries from the FOCUS data spreadsheet recovered.')
 focus_nums_set=set()
 for s in focus_nums.keys():
    focus_nums_set.add(s)
 print("These IDs were not found:", sra_ids.difference(focus_nums_set))
 return(focus_nums)

# 3.
#for each SRA run ID in the dictionary from part 2.:
#get the row of info from the partie data and store in the same dictionary from part 2, appending it to the previous list
def get_partie(focus_nums, partie_runs):
 print("Parsing partie output...")
 sra_ids = set()
 partie_nums={}
 with open(partie_runs,'r') as partie_file:
   for sra_id in focus_nums:
      sra_ids.add(sra_id)
      partie_file.seek(0)
      for line in partie_file:
       if sra_id in line: 
        
        line=line.strip('\n').split('\t')
        partie_id=line
        entries = line[1:6]
        partie_nums[sra_id]=entries

    
 print(len(partie_nums.keys()), ' out of', len(sra_ids),' matching entries from the PARTIE data spreadsheet recovered.')
 partie_nums_set=set()
 for s in partie_nums.keys():
  partie_nums_set.add(s)
 print("These IDs were not found:", sra_ids.difference(partie_nums_set))
 return(partie_nums)

def format_data(partie_nums, focus_nums, run_ids, outf):
 print("Creating output file...")
 sources=set()
 sra_ids_list=list(partie_nums.keys())
 #parse the dictionary to get the full set of sources and percentages
 for f in focus_nums:
    for pair in focus_nums[f].strip('{').strip('}').split(','):
        source= pair.split(':')[0].strip(' ').strip("'")
        sources.add(source)

    #convert the set of sources and the set of sra ids into lists
    #append the partie labels and an environment field  to the header line
 add_labels=['percent_unique_kmer','percent_16S','percent_phage','percent_Prokaryote','percent_human','source']
 header_line=[s for s in sources]
 for p in add_labels:
     header_line.append(p)

    #append the partie numbers to each array
    #create list of lists with the first entry of the list initialized to each SRA-ID
    #and a zero entry for each of the environments
    #so there should be n lists of length e, where n total number of genomes and e is total number of environments
 feature_arrays=[[0] *len(header_line) for i in range(len(sra_ids_list))]

#fill in the percentages from the focus output into the arrays
 for i in range(len(sra_ids_list)):
    entry=focus_nums[sra_ids_list[i]].strip('{').strip('}').split(',')
    for j in range(len(header_line)-6):
        for e in entry:
            pairs=e.split(':')
            source=pairs[0].strip('{').strip(' ').strip("'")
            percentage=float(pairs[1].strip(' ').strip("'"))
            
            if header_line[j] == source:
                    feature_arrays[i][j]=percentage
 for i in range(len(sra_ids_list)):
  entry=partie_nums[sra_ids_list[i]]
  for e in range(0,5):
    feature_arrays[i][len(header_line)-6+e]=entry[e]
 for i in range(len(sra_ids_list)):
     sra=sra_ids_list[i]
     for env in run_ids.keys():
      for s in run_ids[env]:
       if s==sra:
        feature_arrays[i][len(header_line)-1]=env
                    
#add sra ids and sources labels to arrays
 i=0
 for f in feature_arrays:
    f.insert(0,sra_ids_list[i])
    i+=1
 label='sra_id'
 header_line.insert(0,label)
 feature_arrays.insert(0,header_line)
 with open(outf, 'w', newline='') as outf:
  wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
  for f in feature_arrays:
   wr.writerow(f)
 print("Finished!")
 return()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine focus and partie output files and generate a single dataset of  file")
    parser.add_argument('-ids', help='a file of ids and their ground truth values, 2 columns with a header row', required=True)
    parser.add_argument('-focus', help='a json file of the summarized output of the FOCUS tool', required=True)
    parser.add_argument('-partie', help='a tsv file of the summarized output of the partie tool', required=True)
    parser.add_argument('-out', help='output filename', required=True)
    parser.add_argument('-runs', help='the number of runs from each category to include, if not enough runs are found this category will be exlcuded', required=True)
    args = parser.parse_args()
    
    run_ids=get_ids(args.ids)
   # num_runs=int(args.runs)
    focus_nums=get_focus(run_ids, int(args.runs), args.focus)
    partie_nums=get_partie(focus_nums,args.partie)
    format_data(partie_nums, focus_nums, run_ids, args.out)
