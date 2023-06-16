#Part A: cleaning the Data of Kunitz domain downloaded form PDB database (Resolution <3, length 49-80, Sequence Identity grouped by 95%)

tail -n +3 kunitz_pdb.csv | grep -v ",," | cut -d "," -f 2,3 |  tr -d \" | tr "," ":" > kunitz_pdb_idlist.txt  

#Part B: cleaning the  Auth Asym IDs to use them as input for PDB-FOLD website.

tail -n +3 pdb_id.csv | grep -v ",," | cut -d "," -f 2,3 |  tr -d \" | tr "," ":" > pdb_id_clean.txt

#multi align on pdbefold by pdb id list
#remove empty line
grep . kunitz_pdb_msa.fasta | awk 
#
grep . kunitz_pdb_msa.fasta | awk '{if (substr($1,1,1)==">") {printf "\n%s ",$1} else {printf "%s",$0}}' 
#
grep . kunitz_pdb_msa.fasta | awk '{if (substr($1,1,1)==">") {printf "\n%s ",$1} else {printf "%s",$0}}' | awk '{print $1;print substr($2,20,58)}'

#some sequences such as for PDB:5nx1:D are shorter than others so we can exclude them
grep . fasta.seq  | awk '{if (substr($1,1,1)==">") {printf "\n%s ", $1} else {printf "%s",$0}}' | awk '{print $1;print toupper(substr($2,20,58))}' | awk '/^>/{p=1} /5nx1/{p=0} p'|tail -n +3 >ali_3d_clean.fasta

##################################################################
#Model generation
# to build hmm from fasta
hmmbuild aln_3d_clean.hmm aln_3d_clean.fasta 

### The positive and negative set can be downloaded from the uniprot
### with the uniprot ID mapping we make a list of redundant which contains all those kunitz identifiers we used to generate the model with. 
#make a list of redundant
zcat <uniprot-compressed_true_download_true_format_list-2023.05.11-14.26.46.16.list.gz > list_redundant.txt
less list_redundant.txt
wc list_redundant.txt


#Cleaning the all kunitz file (the positive set)and counting
grep ">" kunitz.fasta |wc
grep "^>" kunitz.fasta |cut -d "|" -f 2 >kunitz.list
wc kunitz.list

# comparing all kunitz domain to the list of redundant and removing the intersection
comm <(sort kunitz.list) <(sort list_redundant.txt) |less
comm -23 <(sort kunitz.list) <(sort list_redundant.txt)| wc

##################################################################  
# geting the kunitz fasta file from uniprot after we removed the redundant, #counting it

zcat <uniprot-compressed > kunitz_clean.fasta
grep ">" kunitz_clean.fasta |wc

##################################################################
# Optimization and Model evauation
### Hmmsearch for positve set and cleaning the result.
hmmsearch --max --noali -o kunitz_clean.search aln_3d.hmm kunitz_clean.fasta
#cleaning the kunitz_clean.search file of hmm build for the positive set
head -n 385 kunitz_clean.search
head -n 385 kunitz_clean.search |less
head -n 385 kunitz_clean.search |tail -n +18
head -n 385 kunitz_clean.search |tail -n +18 |less
head -n 385 kunitz_clean.search |tail -n +18 |wc
head -n 385 kunitz_clean.search |tail -n +18 > kunitz_clean.out
# getting the hmm for the negative set
hmmsearch --max --noali -o nonkunitz.search aln_3d.hmm nonkunitz.fasta
# cleaning the nonkunitz.search file of hmm build for the negative set
head -n 57 nonkunitz.search|less
head -n 57 nonkunitz.search|tail -n +18 |less
head -n 57 nonkunitz.search|tail -n +18 |grep -v inclusion |less
head -n 57 nonkunitz.search|tail -n +18 |grep -v inclusion |wc

## Geting the Identifiers and e-value out of the nonkunitz and kunitz file result and assign 0 and 1 to them respectively
awk '{split($9,a,"|"); print a[2],$4,0}' nonkunitz.out >nonkunitz.class

awk '{split($9,a,"|"); print a[2],$4,1}' kunitz_clean.out >kunitz_clean.class

#recover that are missing
grep ">" nonkunitz.fasta|cut -d "|" -f 2 |less
grep ">" nonkunitz.fasta|cut -d "|" -f 2 |wc
grep ">" nonkunitz.fasta|cut -d "|" -f 2 |sort >nonkunitz.list

cut -d " " -f 1 nonkunitz.class |sort >list_hits.txt # number of hits in nonkunitz

comm -23 nonkunitz.list list_hits.txt |awk '{print $0,100,0}' >>nonkunitz.class # missing values

#### two fold cross validation
#make random sets of positive and negative
sort -R nonkunitz.class >nonkunitz.random
sort -R kunitz_clean.class >kunitz_clean.random
head -n 284414 nonkunitz.random >set_1.txt
head -n 184 kunitz_clean.random >>set_1.txt
tail -n +185 kunitz_clean.random >>set_2.txt
tail -n +185 kunitz_clean.random >set_2.txt

# Evaluation of the performance
python3 performance.py set_1.txt # set1 as training set
python3 performance.py set_2.txt # set2 as training set
python3 performance.py (cat set_1.txt  set_2.txt) #evaluate model on both sets