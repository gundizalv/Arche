#!/bin/bash


start_time="$(date -u +%s)"

#Default parameters
#---------------

# Main directory
source ~/.arche.conf
#DIR=/opt/arche_1.0.1


# Procedure
mode=standard
kingdom=bact
kingdom_2="bacteria"
# Name of the created files
input_name=ARCHE
random_append=$(shuf -i1-100000 -n1)
# Search for bacterial rRNAs
barrnap_mode=bac
# Search for bacterial tRNAs
tRNA_mode=-B
# Genome type option for GeneMark2-S
gms2_mode=bacteria
# Uniprot database
uniprot_db=bac
uni_db_kingdom=
# Number of CPUs to use
cpus=1
# Similarity e-value cut-off
evalue=1e-08
# Amount of RAM to use in Gb
memory=2
# Structural annotator. Don't change this!
annotator=genemark
db_search_alig=blastp
# Path of the output directory.
dir_path=
current_dir=$(realpath .)
# Minimum coverage on query protein
cov=70
# Bypass RNA prediction
bypass_2="no"

# Console format (colour and bold)
bold=$(tput bold)
cyan=$(tput setaf 6)
yellow=$(tput setaf 3)
red=$(tput setaf 1)
green=$(tput setaf 2)
reset=$(tput sgr0)

#Verbose options
verbose_on=
quiet_on=--quiet
quiet_on_prodigal=-q
verbose_2="off"


#-------------------------------------------------------------------------------------------------------------------------------------------------
#Functions.
#-------------------------------------------------------------------------------------------------------------------------------------------------
check_success() {
	if [[ ! $? -eq 0 ]]; then 
		echo "
ERROR: Command failed — aborting
		"
		exit 1
	fi

}

validate_numeric() {
	if [[ ! "$2" =~ ^[0-9]+([.][0-9]+)?$ && \
		  ! "$2" =~ ^[0-9]+([.][0-9]+)?e-[0-9]+$ ]]; then
		echo "
ERROR: Option '$1' requires a numeric value. You provided '$2'
"
		exit 1
	fi

}

usage() {
echo "
Arche 1.0.0
---------
Written by Daniel Alonso <gundizalvus16@hotmail.com>
Last modified May 7, 2025
Description:
        Optimized annotator for prokariotic contigs.

Usage:
	arche --download
	arche --install
        arche --help
	arche [options] <contigs.fasta>

Options:
	-h, --help           This help.
	-d, --download	     Download and extract Arche's database directory into the current directory.
	-i, --install	     Set up the Arche's database directory location, and install databases.
	-n, --name-files     Name of the files to be created in the output directory, including the
			     directory itself. Default 'arche'.
	-o, --output	     Provide the full path to the directory where the output directory will be
			     created. E.g. /home/user/ . Default current.
	-k, --kingdom        Source of the contigs. Use 'arch' for archaeal genomes or 'meta' for metagenomes
			     (default is for bacterial genomes).
	-m, --mode           Gives priority to Orthology (KO, eggNOG) or Enzyme Comission designed databases
			     during the annotation. Use 'kegg' for KO-->eggNOG-->E.C., 'eggnog' for
			     eggNOG-->KO-->E.C., or 'ec' for E.C.-->KO-->eggNOG (default will use a shorter
			     swiss-prot KO·eggNOG·E.C. designed database with no priority).
	-a, --alignment      Select the algorithm to use during the protein alignment step: 'diamond'
			     (accelerated blastp) or 'ssearch' (Smith-Waterman). Default 'blastp'.
	-t, --threads        Number of threads to use (default '1').
	-r, --memory         Amount of RAM to use in GB (default '2').
	-e, --evalue         Similarity e-value cut-off (default '1e-08').
	-q, --query-cov      Minimum coverage on query protein (default '70').
	-b, --bypass         Use 'yes' to bypass the RNA gene prediction.
	-v, --verbose	     Use 'yes' to turn on the verbose mode.
	"
}

download_db(){
	rm arche_*.tar 2>/dev/null
	gdown --fuzzy https://drive.google.com/file/d/1x9caXGPpYXCHUoodOdnuJI0tCDe9qtGG/view?usp=sharing
	check_success
	echo "
Extracting..."
	tar xvf arche_*.tar
	check_success
	echo "
Done.
Move Arche's database directory to the desired location and run arche.sh --install"
}

install_dir () {
	echo "
Enter the full path to the Arche's database directory.
"
	read -e db_path

	if [[ ! -d "$db_path" ]]; then
		echo "ERROR: '$db_path' is not a valid directory"
		exit 1
	fi

	echo "DIR=$db_path" > ~/.arche.conf
	echo "Database path saved to .arche.conf"
	source ~/.arche.conf

	#main_dir=$(realpath arche.sh | sed 's!/bin/arche.sh!!')
	#sed -i -z "s|DIR=[^\n]*|DIR=$main_dir|" arche.sh
}

install_db () {
echo "
${bold}${cyan}Unpacking databases and map files...${reset}
"
7z e "$DIR"/db/nfam/NF.hmm.7z -o"$DIR"/db/nfam/
check_success
7z e "$DIR"/db/ssearch/EC_db.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/EC_db_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/EC_db_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG-EC_db.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG-EC_db_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG-EC_db_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/eggNOG-EC_db.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/eggNOG-EC_db_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/eggNOG-EC_db_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/sprot.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/sprot_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/sprot_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG_db.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG_db_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG_db_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG-eggNOG_db.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG-eggNOG_db_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/KEGG-eggNOG_db_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/eggNOG_db.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/eggNOG_db_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/eggNOG_db_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/EC-KEGG-eggNOG_db.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_MET.fa.7z -o"$DIR"/db/ssearch/
7z e "$DIR"/db/kofam/KO_prokaryote.hmm.7z -o"$DIR"/db/kofam/
7z e "$DIR"/db/tigrfam/TIGRFAM.hmm.7z -o"$DIR"/db/tigrfam/
7z e "$DIR"/db/pfam/pfam-A.hmm.7z -o"$DIR"/db/pfam/
7z e "$DIR"/mpp/NFAM_EC3.tsv.7z -o"$DIR"/mpp/
7z e "$DIR"/mpp/uniprot_bac.7z -o"$DIR"/mpp/
7z e "$DIR"/mpp/uniprot_arch.7z -o"$DIR"/mpp/
7z e "$DIR"/mpp/uniprot_met.7z -o"$DIR"/mpp/
7z e "$DIR"/mpp/KO_match_table6.tsv.7z -o"$DIR"/mpp/
7z e "$DIR"/mpp/TIGRFAM_EC3.tsv.7z -o"$DIR"/mpp/
7z e "$DIR"/mpp/PFAM_EC7.tsv.7z -o"$DIR"/mpp/

echo "
${bold}${cyan}Setting up databases...${reset}
"
diamond makedb --in "$DIR"/db/ssearch/sprot.fa -d "$DIR"/db/diamond/sprot
check_success
diamond makedb --in "$DIR"/db/ssearch/sprot_ARCH.fa -d "$DIR"/db/diamond/sprot_ARCH
diamond makedb --in "$DIR"/db/ssearch/sprot_MET.fa -d "$DIR"/db/diamond/sprot_MET
diamond makedb --in "$DIR"/db/ssearch/EC_db.fa -d "$DIR"/db/diamond/EC_db
diamond makedb --in "$DIR"/db/ssearch/EC_db_ARCH.fa -d "$DIR"/db/diamond/EC_db_ARCH
diamond makedb --in "$DIR"/db/ssearch/EC_db_MET.fa -d "$DIR"/db/diamond/EC_db_MET
diamond makedb --in "$DIR"/db/ssearch/KEGG_db.fa -d "$DIR"/db/diamond/KEGG_db
diamond makedb --in "$DIR"/db/ssearch/KEGG_db_ARCH.fa -d "$DIR"/db/diamond/KEGG_db_ARCH
diamond makedb --in "$DIR"/db/ssearch/KEGG_db_MET.fa -d "$DIR"/db/diamond/KEGG_db_MET
diamond makedb --in "$DIR"/db/ssearch/eggNOG_db.fa -d "$DIR"/db/diamond/eggNOG_db
diamond makedb --in "$DIR"/db/ssearch/eggNOG_db_ARCH.fa -d "$DIR"/db/diamond/eggNOG_db_ARCH
diamond makedb --in "$DIR"/db/ssearch/eggNOG_db_MET.fa -d "$DIR"/db/diamond/eggNOG_db_MET
diamond makedb --in "$DIR"/db/ssearch/KEGG-eggNOG_db.fa -d "$DIR"/db/diamond/KEGG-eggNOG_db
diamond makedb --in "$DIR"/db/ssearch/KEGG-eggNOG_db_ARCH.fa -d "$DIR"/db/diamond/KEGG-eggNOG_db_ARCH
diamond makedb --in "$DIR"/db/ssearch/KEGG-eggNOG_db_MET.fa -d "$DIR"/db/diamond/KEGG-eggNOG_db_MET
diamond makedb --in "$DIR"/db/ssearch/KEGG-EC_db.fa -d "$DIR"/db/diamond/KEGG-EC_db
diamond makedb --in "$DIR"/db/ssearch/KEGG-EC_db_ARCH.fa -d "$DIR"/db/diamond/KEGG-EC_db_ARCH
diamond makedb --in "$DIR"/db/ssearch/KEGG-EC_db_MET.fa -d "$DIR"/db/diamond/KEGG-EC_db_MET
diamond makedb --in "$DIR"/db/ssearch/eggNOG-EC_db.fa -d "$DIR"/db/diamond/eggNOG-EC_db
diamond makedb --in "$DIR"/db/ssearch/eggNOG-EC_db_ARCH.fa -d "$DIR"/db/diamond/eggNOG-EC_db_ARCH
diamond makedb --in "$DIR"/db/ssearch/eggNOG-EC_db_MET.fa -d "$DIR"/db/diamond/eggNOG-EC_db_MET
diamond makedb --in "$DIR"/db/ssearch/EC-KEGG-eggNOG_db.fa -d "$DIR"/db/diamond/EC-KEGG-eggNOG_db
diamond makedb --in "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa -d "$DIR"/db/diamond/EC-KEGG-eggNOG_db_ARCH
diamond makedb --in "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_MET.fa -d "$DIR"/db/diamond/EC-KEGG-eggNOG_db_MET

makeblastdb -in "$DIR"/db/ssearch/sprot.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/sprot
check_success
makeblastdb -in "$DIR"/db/ssearch/sprot_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/sprot_ARCH
makeblastdb -in "$DIR"/db/ssearch/sprot_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/sprot_MET
makeblastdb -in "$DIR"/db/ssearch/EC_db.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/EC_db
makeblastdb -in "$DIR"/db/ssearch/EC_db_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/EC_db_ARCH
makeblastdb -in "$DIR"/db/ssearch/EC_db_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/EC_db_MET
makeblastdb -in "$DIR"/db/ssearch/KEGG_db.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG_db
makeblastdb -in "$DIR"/db/ssearch/KEGG_db_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG_db_ARCH
makeblastdb -in "$DIR"/db/ssearch/KEGG_db_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG_db_MET
makeblastdb -in "$DIR"/db/ssearch/eggNOG_db.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/eggNOG_db
makeblastdb -in "$DIR"/db/ssearch/eggNOG_db_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/eggNOG_db_ARCH
makeblastdb -in "$DIR"/db/ssearch/eggNOG_db_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/eggNOG_db_MET
makeblastdb -in "$DIR"/db/ssearch/KEGG-eggNOG_db.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG-eggNOG_db
makeblastdb -in "$DIR"/db/ssearch/KEGG-eggNOG_db_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG-eggNOG_db_ARCH
makeblastdb -in "$DIR"/db/ssearch/KEGG-eggNOG_db_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG-eggNOG_db_MET
makeblastdb -in "$DIR"/db/ssearch/KEGG-EC_db.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG-EC_db
makeblastdb -in "$DIR"/db/ssearch/KEGG-EC_db_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG-EC_db_ARCH
makeblastdb -in "$DIR"/db/ssearch/KEGG-EC_db_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/KEGG-EC_db_MET
makeblastdb -in "$DIR"/db/ssearch/eggNOG-EC_db.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/eggNOG-EC_db
makeblastdb -in "$DIR"/db/ssearch/eggNOG-EC_db_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/eggNOG-EC_db_ARCH
makeblastdb -in "$DIR"/db/ssearch/eggNOG-EC_db_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/eggNOG-EC_db_MET
makeblastdb -in "$DIR"/db/ssearch/EC-KEGG-eggNOG_db.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/EC-KEGG-eggNOG_db
makeblastdb -in "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/EC-KEGG-eggNOG_db_ARCH
makeblastdb -in "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_MET.fa -dbtype prot -parse_seqids -out "$DIR"/db/blastp/EC-KEGG-eggNOG_db_MET

map_db "$DIR"/db/ssearch/sprot.fa
check_success
map_db "$DIR"/db/ssearch/sprot_ARCH.fa
map_db "$DIR"/db/ssearch/sprot_MET.fa
map_db "$DIR"/db/ssearch/EC_db.fa
map_db "$DIR"/db/ssearch/EC_db_ARCH.fa
map_db "$DIR"/db/ssearch/EC_db_MET.fa
map_db "$DIR"/db/ssearch/KEGG_db.fa
map_db "$DIR"/db/ssearch/KEGG_db_ARCH.fa
map_db "$DIR"/db/ssearch/KEGG_db_MET.fa
map_db "$DIR"/db/ssearch/eggNOG_db.fa
map_db "$DIR"/db/ssearch/eggNOG_db_ARCH.fa
map_db "$DIR"/db/ssearch/eggNOG_db_MET.fa
map_db "$DIR"/db/ssearch/KEGG-eggNOG_db.fa
map_db "$DIR"/db/ssearch/KEGG-eggNOG_db_ARCH.fa
map_db "$DIR"/db/ssearch/KEGG-eggNOG_db_MET.fa
map_db "$DIR"/db/ssearch/KEGG-EC_db.fa
map_db "$DIR"/db/ssearch/KEGG-EC_db_ARCH.fa
map_db "$DIR"/db/ssearch/KEGG-EC_db_MET.fa
map_db "$DIR"/db/ssearch/eggNOG-EC_db.fa
map_db "$DIR"/db/ssearch/eggNOG-EC_db_ARCH.fa
map_db "$DIR"/db/ssearch/eggNOG-EC_db_MET.fa
map_db "$DIR"/db/ssearch/EC-KEGG-eggNOG_db.fa
map_db "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa
map_db "$DIR"/db/ssearch/EC-KEGG-eggNOG_db_MET.fa

hmmpress "$DIR"/db/kofam/KO_prokaryote.hmm
check_success
hmmpress "$DIR"/db/pfam/pfam-A.hmm
hmmpress "$DIR"/db/tigrfam/TIGRFAM.hmm
hmmpress "$DIR"/db/nfam/NF.hmm

rm "$DIR"/db/nfam/NF.hmm.7z
rm "$DIR"/db/ssearch/*.7z
rm "$DIR"/db/kofam/KO_prokaryote.hmm.7z
rm "$DIR"/db/tigrfam/TIGRFAM.hmm.7z
rm "$DIR"/db/pfam/pfam-A.hmm.7z
rm "$DIR"/mpp/*.7z
}

heurist_search () {
	case "$db_search_alig" in
	blastp)
		index=
		useless_lines=#
		cut_col=2
		cut_ID=3
		cut_table=5
		select_db
		blastp -db "$db_path" -query "$input_fasta" -evalue "$evalue" -out "$database"_out -max_target_seqs 1 -num_threads "$cpus" -outfmt "7 qseqid sseqid qcovs"
		check_success
		match=$(grep -v "#" "$database"_out | awk -v cov2=$cov '{if ($3 >= cov2) print $1}' | sort -V | uniq -u)
		[ -n "$match" ] && heurist_func_annot
		report
		;;
	diamond)
		let block_size=$memory/6
		index=.dmnd
		useless_lines=*
		cut_col=2
		cut_ID=2
		cut_table=3
		select_db
		diamond blastp $quiet_on $verbose_on --db "$db_path" --query "$input_fasta" --sensitive --outfmt 6 --max-target-seqs 1 --evalue "$evalue" --query-cover "$cov" --unal 1 --block-size "$block_size" --out "$database"_out
		check_success
		match=$(cat "$database"_out | grep -v "*" | cut -f1 | sort -V | uniq -u)
		[ -n "$match" ] && heurist_func_annot
		report
		;;
	ssearch)
		index=.fa
		useless_lines=#
		cut_col=3
		cut_ID=2
		cut_table=3
		select_db	
		ssearch36 -m8CBl -b 1 -E "$evalue" "$input_fasta" "$db_path" > "$database"_out
		check_success
		match=$(grep -v "#" "$database"_out | awk -v cov2=$cov '{if ((($6/$2)*100 >= cov2)) print $1"\t"$3}' | cut -f1 | sort -V | uniq -u)	
		[ -n "$match" ] && heurist_func_annot
		report
		;;	
	esac
} 

select_db () {
	[ "$db_heurist" = dbST ] && db_path="$DIR"/db/"$db_search_alig"/sprot"$uni_db_kingdom""$index"
	[ "$db_heurist" = db00 ] && db_path="$DIR"/db/"$db_search_alig"/EC-KEGG-eggNOG_db"$uni_db_kingdom""$index" #EC-KEGG-eggNOG
	[ "$db_heurist" = db01 ] && db_path="$DIR"/db/"$db_search_alig"/KEGG-eggNOG_db"$uni_db_kingdom""$index" #KEGG-eggNOG
	[ "$db_heurist" = db02 ] && db_path="$DIR"/db/"$db_search_alig"/KEGG-EC_db"$uni_db_kingdom""$index" #KEGG-EC
	[ "$db_heurist" = db03 ] && db_path="$DIR"/db/"$db_search_alig"/KEGG_db"$uni_db_kingdom""$index" #KEGG
	[ "$db_heurist" = db04 ] && db_path="$DIR"/db/"$db_search_alig"/eggNOG-EC_db"$uni_db_kingdom""$index" #eggNOG-EC.tsv
	[ "$db_heurist" = db05 ] && db_path="$DIR"/db/"$db_search_alig"/eggNOG_db"$uni_db_kingdom""$index" #eggNOG.tsv
	[ "$db_heurist" = db06 ] && db_path="$DIR"/db/"$db_search_alig"/EC_db"$uni_db_kingdom""$index" #EC.tsv
}

heurist_func_annot () {
	cat "$database"_out | grep -v "$useless_lines" | cut -f1,"$cut_col" | sed -E 's/\t/|/g' | sed 's/$/|/g' > MAP_LIST
	cut -f"$cut_ID" -d "|" MAP_LIST  > new_LIST
	faSomeRecords "$DIR"/mpp/uniprot_"$uniprot_db" new_LIST OUT 
	grep ">" OUT | sed 's/>//g' > BANK
	cat new_LIST | while read -r line ; do grep -w -m 1 "$line" BANK >> post_TABLE || printf "\n" >> post_TABLE ; done
	paste MAP_LIST post_TABLE | cut -f1,"$cut_table"-12 -d "|" > "$database"_TABLE
	rm OUT BANK MAP_LIST new_LIST post_TABLE
}

report () {
	if [[ -n "$match" ]]; then
		m_count=$(echo "$match" | wc -l)
		echo -n "${bold}${cyan}$m_count match(es) reported${reset}"
		echo "$match" > match_"$database"_list.txt
		faSomeRecords -exclude "$input_fasta" match_"$database"_list.txt "$input_name"_"$database"_non_match.faa
		rm match_"$database"_list.txt
	else
		m_count=0
		echo -n "${bold}${cyan}0 match(es) reported${reset}"
		cp "$input_fasta" "$input_name"_"$database"_non_match.faa
	fi
}

cite_blastp () {
	echo "
${bold}Querying Uniprot with BLASTp+ (heuristic algorithm)${reset} — The UniProt Consortium. UniProt: the universal protein knowledgebase in 2021 (2021). Nucleic Acids Res. 49:D1; Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421.
"
}

cite_diamond () { 
	echo "
	
${bold}Queriying Uniprot with DIAMOND (fast heuristic algorithm)${reset} — The UniProt Consortium. UniProt: the universal protein knowledgebase in 2021 (2021). Nucleic Acids Res. 49:D1; Buchfink B, Xie C, Huson DH, 'Fast and sensitive protein alignment using DIAMOND', Nature Methods 12, 59-60 (2015). doi:10.1038/nmeth.3176
"
}

cite_ssearch () {
echo "
${bold}Queriying Uniprot with SSEARCH (Smith-Waterman algorithm)${reset} — The UniProt Consortium. UniProt: the universal protein knowledgebase in 2021 (2021). Nucleic Acids Res. 49:D1; W. R. Pearson. Effective protein sequence comparison (1996). Methods Enzymol., 266:227–258
"
}

kofam_search () {
	echo "
	
${bold}Querying KofamKOALA (prokaryotic) with HMMER3${reset} — S. R. Eddy. (2011)  Accelerated profile HMM searches. PLoS Comp. Biol., 7:e1002195.; Takuya Aramaki, Romain Blanc-Mathieu, Hisashi Endo, Koichi Ohkubo, Minoru Kanehisa, Susumu Goto, Hiroyuki Ogata. (2020) KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 36:2251–2252
"
	input_fasta="$input_name"_"$database"_non_match.faa
	if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_kofam_out "$DIR"/db/kofam/KO_prokaryote.hmm "$input_fasta" #HERE ADDED '-E "$evalue" --domE "$evalue" -Z 10939 --domZ 10939'. Then removed -E "$evalue" --domE "$evalue" -Z 10939 --domZ 10939
	else
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_kofam_out "$DIR"/db/kofam/KO_prokaryote.hmm "$input_fasta" > log_hmm_kofam	#HERE ADDED '-E "$evalue" --domE "$evalue" -Z 10939 --domZ 10939'. Then removed -E "$evalue" --domE "$evalue" -Z 10939 --domZ 10939
	fi
	check_success
	hmmscan_out=$(<hmmscan_kofam_out)
	database=kofam
	map_file=KO_match_table6.tsv
	cut_column=1
	hmm_func_annot
	match=$(echo "$hmmscan_out" | list_hmm_hits)
	report ; echo "${bold}${cyan} with KEGG only associated number(s) (from HMMs).${reset}"
	echo "$m_count KEGG (from HMMs)" >> arche_report
}

hmm_func_annot () {
	echo "$hmmscan_out" | sed '/^#/d' | awk '!x[$4]++' | awk -v cov3=$cov -v column=$cut_column '{if ((($17-$16)/$6)*100 >= cov3) print $column"\t"$4}' | sed -E 's/\t/|/g' |  sed 's/$/|/g' | sort  > pre_TABLE
	map_LIST=$(cut -f1 -d "|" pre_TABLE | cut -f1 -d ".")
	echo "$map_LIST" | while read -r line ; do grep -w -m 1 "$line" "$DIR"/mpp/"$map_file" >> post_TABLE || printf "\n" >> post_TABLE ; done
	paste pre_TABLE post_TABLE  | sed -E 's/\t//g' | cut -d "|" -f2-10 > "$database"_TABLE
	rm pre_TABLE post_TABLE
}

list_hmm_hits () {
	sed '/^#/d' | awk '!x[$4]++' | awk -v cov3=$cov '{if ((($17-$16)/$6)*100 >= cov3) print $1"\t"$4}' | cut -f2 | sort -V | uniq -u
}



#-------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------


#Options
#------

options=$(getopt -l "install,download,kingdom:,bypass:,mode;,alignment:,name-files:,output:,threads:,evalue:,memory:,query-cov:,help,verbose:" -o "idk:b:m:a:n:o:t:e:r:q:hv:" -- "$@")
[ $? -eq 0 ] || {
	echo "Incorrect options provided."
    exit 1
}
eval set -- "$options"
while true ; do
	case $1 in
	-i|--install)
		install_dir
		install_db
		exit 0
		;;
	-d|--download)
		download_db
		exit 0
		;;
	-k|--kingdom)
		shift
		[[ ! "$1" =~ ^(arch|meta)$ ]] && {
            echo "Incorrect option provided for '-k'/'--kingdom'. Choose 'arch' for archaea or 'meta' for metagenome. Remove the option for bacterial genomes."
            exit 1
        }
		 kingdom=$1
		;;
	-b|--bypass)
		shift; 
		if [[ "$1" = yes ]]; then
			bypass="true"
			bypass_2="yes"
		else
			echo "Wrong argument to '--bypass' option. Use 'yes' to bypass the RNA gene prediction"
    	exit 1	
		fi	
		;;
	-m|--mode)
		shift
		[[ ! "$1" =~ ^(kegg|eggnog|ec)$ ]] && {
            echo "Incorrect option provided for '-m'/'--mode'. Choose 'kegg', 'eggnog' or 'ec'."
            exit 1
        }
		 mode=$1
		;;
	-a|--alignment)
		shift
		[[ ! "$1" =~ ^(diamond|ssearch)$ ]] && {
            echo "Incorrect option provided for '-a'/'--engine'. Choose 'diamond' or 'ssearch'."
            exit 1
        }
		 db_search_alig=$1
		;;
	-n|--name-files)
		shift
		input_name=$1
		;;
	-o|--output)
		shift
		if [[ "$1" == */ ]]; then
			dir_path=$1
		else
			dir_path="$1"/
		fi
		[ -d "$dir_path" ] || {
		echo "
ERROR: wrong output path provided. Make sure it ends with a slash (/) and it doesn't includes the name of the output directory. Check the help for an example.
		"
    	exit 1
		}
		;;
	-t|--threads)
		shift
		cpus=$1
		validate_numeric "--threads" "$cpus"
		;;
	-e|--evalue)
		shift
		evalue=$1
		validate_numeric "--evalue" "$evalue"
		;;
	-r|--memory)
		shift
		memory=$1
		validate_numeric "--memory" "$memory"
		;;
	-q|--query-cov)
		shift
		cov=$1
		validate_numeric "--query-cov" "$cov"
		;;
	-h|--help)
		usage
		exit 0
		;;
	-v|--verbose)
		shift; 
		if [[ "$1" = yes ]]; then
			quiet_on=
			quiet_on_prodigal=
			verbose_on=-v
			verbose_2="on"
			set -x
		else
			echo "Wrong argument to '--verbose' option. Use 'yes' to activate verbose mode"
    	exit 1	
		fi	
		
		;;
	--)
		shift
	    break;;
	esac
	shift
done

# -------- GeneMarkS-2 validation --------
if ! command -v gms2.pl >/dev/null 2>&1; then
	echo "
ERROR: GeneMarkS-2 (gms2.pl) is not installed, not in your PATH, or its license has expired.

Please install it from http://exon.gatech.edu/GeneMark/license_download.cgi
"
	exit 1
fi

# -------- Arche's database directory validation --------

if [[ -z "$DIR" || ! -d "$DIR" ]]; then
    echo "ERROR: The variable DIR is not set or does not point to a valid directory."
    echo "Please run: arche.sh --install"
    exit 1
fi

# -------- Input FASTA file validation --------
# Ensures a valid file with .fa, .fna, .fasta or .ffn extension is provided as last argument

if [[ -z "${!#}" ||
	! -f "${!#}" ||
	! "${!#}" =~ [.](fa|fna|fasta|ffn)$ ]]; then
	echo "ERROR: You must provide a valid FASTA file (.fa, .fna, .fasta, or .ffn)."
	echo "Example: arche.sh -n sample -t 4 sample_genome.fna"
	exit 1
fi



# In case of archaeal genome, or metagenome
#------------------------------------

if [[ "$kingdom" = arch ]]; then	
	 barrnap_mode=arc
	 tRNA_mode=-A
	 gms2_mode=archaea
	uniprot_db=arch
	uni_db_kingdom=_ARCH
	kingdom_2="archaea"
elif [[ "$kingdom" = meta ]]; then	
	 annotator=metaprodigal
	 db_search_alig=diamond
	uniprot_db=met
	uni_db_kingdom=_MET
	kingdom_2="metagenome"
	bypass_2="yes"
fi

# Parameters of the run
#--------------------

echo	"
${bold}${cyan}Runing Arche with the following parameters:${reset}

${bold}input = "${!#}"
kingdom = "$kingdom_2"
mode = "$mode"
search algorithms = "$db_search_alig" and HMMER3
cpus = "$cpus"
e-value = "$evalue"
memory = "$memory" GB
bypass RNA prediction = "$bypass_2"
gene prediction = "$annotator"
query coverage = "$cov" %
verbose mode = "$verbose_2"${reset}
"


# Create output directory
# -------------------

directory_name=$(echo ${input_name^^} | cut -c1-8)
if [[ -d "$dir_path""$directory_name"_"$random_append" ]]; then
	echo "
${bold}${yellow}The directory "$directory_name"_"$random_append" exists! Remove or rename it.${reset}
"
	exit 1
else
	echo "
${bold}Creating output directory "$directory_name"_"$random_append"${reset}"
	mkdir  "$dir_path""$directory_name"_"$random_append"
	cp "${!#}" "$dir_path""$directory_name"_"$random_append"/
	cd "$dir_path""$directory_name"_"$random_append"
fi
	
echo "Runing Arche with the following parameters:" >> arche_report
echo "input = ${!#}" >> arche_report
echo "kingdom = $kingdom_2" >> arche_report
echo "mode = $mode" >> arche_report
echo "cpus = $cpus" >> arche_report
echo "e-value = $evalue" >> arche_report
echo "memory = $memory GB" >> arche_report
echo "bypass RNA prediction = $bypass_2" >> arche_report
echo "gene prediction = $annotator" >> arche_report
echo "query coverage = $cov %" >> arche_report
echo "verbose mode = $verbose_2" >> arche_report



# Check contig name for further "sort -V"
#----------------------------------
if [[ "$kingdom" != meta ]]; then
	contig_file="${!#}"
	contig_name_start=$(grep ">" "$contig_file" | head -n 1 | sed 's/>//' | cut -c1);
	if [[ "$contig_name_start" != [A-Za-z] ]]; then
		echo "
${bold}${yellow}Contig headers start with non-letter characters... creating a new contig file ("$contig_file"_2). ${reset}"
		sed 's/>/>c_/g' "$contig_file" > "$contig_file"_2
		contig_file="$contig_file"_2
	fi
fi




if [[ "$bypass" != "true" ]]; then

# Ribosomal RNA gene prediction.
#--------------------------

	[ $annotator = genemark ] && {	
		echo "

${bold}Predicting ribosomal RNA with Barrnap${reset} — Seemann T. barrnap 0.9: rapid ribosomal RNA prediction https://github.com/tseemann/barrnap
"
		barrnap $quiet_on --kingdom  "$barrnap_mode" --outseq rRNA.fna --threads "$cpus" "$contig_file" > rRNA.tsv
		check_success
		barrnap_exit=$(grep ">" rRNA.fna | wc -l)
		echo "
${bold}${cyan}$barrnap_exit rRNAs gene(s) reported.${reset}"
}	

echo "" >> arche_report
echo "Results:" >> arche_report
echo "$barrnap_exit rRNAs" >> arche_report


# tRNA gene prediction.
#------------------

	[ $annotator = genemark ] && {
		echo "
${bold}Predicting transfer RNA with tRNAscan-SE${reset} — Chan, P.P. and Lowe, T. M. (2019) tRNAscan-SE: Searching for tRNA Genes in Genomic Sequences. Methods Mol Biol. 1962:1-14.
"
		tRNAscan-SE $tRNA_mode $quiet_on -o tRNA.tsv -a tRNA.fna "$contig_file"
		check_success
		tRNAscan_exit=$(grep ">" tRNA.fna | wc -l)
		echo "
${bold}${cyan}$tRNAscan_exit tRNAs gene(s) reported.${reset}"
		echo "$tRNAscan_exit tRNAs" >> arche_report
}

fi



# Structural annotation.
#------------------

if [[ $annotator = genemark ]]; then 	
	echo "
${bold}Structural annotation with GeneMarkS-2${reset} — Lomsadze A, Gemayel K, Tang S, Borodovsky M (2018). Modeling leaderless transcription and atypical genes results in more accurate gene prediction in prokaryotes. Genome Res, 29(7), pp 1079-1089
"	
	gms2.pl --gcode 11 --seq  "$contig_file" --fnn "$input_name"_struc_annot.fna --faa "$input_name"_struc_annot.faa  --genome-type "$gms2_mode" --output "$input_name"_GeneMark.lst
	check_success
	features_tot=$(cat "$input_name"_struc_annot.faa | grep ">" | wc -l)
	atypical_tot=$(cat "$input_name"_struc_annot.faa | grep "atypical" | wc -l)
	let native=$features_tot-$atypical_tot
	echo "${bold}${cyan}$features_tot feature(s) reported; $native were native and $atypical_tot atypical.${reset}"
	echo "$features_tot putative protein(s) ($native native and $atypical_tot atypical)" >> arche_report
else
	echo "

${bold}Structural annotation with MetaProdigal ${reset} — Hyatt,D., LoCascio,P.F., Hauser,L.J. and Uberbacher,E.C. (2012) Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics, 28, 2223–2230
" 
	prodigal -g 11 -a "$input_name"_struc_annot.faa -d "$input_name"_struc_annot.fna -f gff -o "$input_name"_struc_annot.gff $quiet_on_prodigal -p meta -i "${!#}" > prodigal.log
	check_success
	features_tot=$(cat "$input_name"_struc_annot.faa | grep ">" | wc -l)
	echo "
${bold}${cyan}$features_tot feature(s) reported.${reset}"
	echo "$features_tot putative protein(s)." >> arche_report
fi



#Maximizing KEGG, eggNOG, E.C. or use SwissProt database only (standard).
#--------------------------------------------------------------

case "$mode" in
standard)
	# Protein database searching with heuristic/SW methods.
	cite_"$db_search_alig"
	input_fasta="$input_name"_struc_annot.faa
	database=heuristic
	db_heurist=dbST
	heurist_search
	echo "$m_count standard (SwissProt) match(es)" >> arche_report
	# Kofam database.
	kofam_search
	;;
	
kegg)	
	# Prioritize target sequences with KEGG associated IDs.
	cite_"$db_search_alig"
	input_fasta="$input_name"_struc_annot.faa
	database=heuristic0
	db_heurist=db00 #EC-KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with KEGG + eggNOG + E.C. associated number(s).${reset}"
	echo "$m_count KEGG + eggNOG + E.C." >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic1
	db_heurist=db01 #KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with KEGG + eggNOG associated number(s).${reset}"
	echo "$m_count KEGG + eggNOG" >> arche_report

	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic2
	db_heurist=db02 #KEGG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with KEGG + E.C. associated number(s).${reset}"
	echo "$m_count KEGG + E.C." >> arche_report	

	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic3
	db_heurist=db03 #KEGG.tsv
	heurist_search ; echo "${bold}${cyan} with KEGG only associated number(s) (from alignments).${reset}"
	echo "$m_count KEGG (from alignments)" >> arche_report	

	# Kofam database.
	kofam_search ; echo ""
	
	# ... continuation.
	cite_"$db_search_alig"
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic4
	db_heurist=db04 #eggNOG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + E.C. associated number(s).${reset}"
	echo "$m_count eggNOG + E.C." >> arche_report

	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic5
	db_heurist=db05
	heurist_search ; echo "${bold}${cyan} with eggNOG only associated number(s).${reset}"
	echo "$m_count eggNOG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic6
	db_heurist=db06 #EC.tsv
	heurist_search ; echo "${bold}${cyan} with E.C. only associated number(s).${reset}"
	echo "$m_count E.C." >> arche_report


	;;
eggnog)
	# Prioritize target sequences with eggNOG associated IDs.
	cite_"$db_search_alig"
	input_fasta="$input_name"_struc_annot.faa
	database=heuristic0
	db_heurist=db00 #EC-KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + KEGG + E.C. associated number(s).${reset}"
	echo "$m_count eggNOG + KEGG + E.C." >> arche_report

	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic1
	db_heurist=db01 #KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + KEGG associated number(s).${reset}"
	echo "$m_count eggNOG + KEGG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic2
	db_heurist=db04 #eggNOG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + E.C. associated number(s).${reset}"
	echo "$m_count eggNOG + E.C." >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic3
	db_heurist=db05 #eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG only associated number(s).${reset}"
	echo "$m_count eggNOG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic4
	db_heurist=db02 #KEGG-EC
	heurist_search ; echo "${bold}${cyan} with KEGG + E.C. associated number(s).${reset}"
	echo "$m_count KEGG + E.C." >> arche_report

	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic5
	db_heurist=db03 #KEGG
	heurist_search ; echo "${bold}${cyan} with KEGG only associated number(s) (from alignments).${reset}"
	echo "$m_count KEGG (from alignments)" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic6
	db_heurist=db06 #EC
	heurist_search ; echo "${bold}${cyan} with E.C. only associated number(s).${reset}"
	echo "$m_count E.C." >> arche_report

	
	# Kofam database.
	kofam_search
	;;
	
ec)
	# Prioritize target sequences with associated E.C. numbers.
	cite_"$db_search_alig"
	input_fasta="$input_name"_struc_annot.faa
	database=heuristic0
	db_heurist=db00 #EC-KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with E.C. + KEGG + eggNOG associated number(s).${reset}"
	echo "$m_count E.C. + KEGG + eggNOG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic1
	db_heurist=db02 #KEGG-EC
	heurist_search ; echo "${bold}${cyan} with E.C. + KEGG associated number(s).${reset}"
	echo "$m_count E.C. + KEGG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic2
	db_heurist=db04 #eggNOG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with E.C. + eggNOG associated number(s).${reset}"
	echo "$m_count E.C. + eggNOG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic3
	db_heurist=db06 #EC.tsv
	heurist_search ; echo "${bold}${cyan} with E.C. only associated number(s).${reset}"
	echo "$m_count E.C." >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic4
	db_heurist=db01 #KEGG-eggNOG
	heurist_search ; echo "${bold}${cyan} with KEGG + eggNOG associated number(s).${reset}"
	echo "$m_count KEGG + eggNOG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic5
	db_heurist=db03 #KEGG
	heurist_search ; echo "${bold}${cyan} with KEGG only associated number(s) (from alignments).${reset}"
	echo "$m_count KEGG" >> arche_report
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic6
	db_heurist=db05 #eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG only associated number(s).${reset}"
	echo "$m_count eggNOG" >> arche_report
	
	# Kofam database.
	kofam_search
	;;
	
esac



# NFAM database.
#--------------

echo "

${bold}Querying NCBI Protein Family Models with HMMER3${reset} — Li W, O'Neill KR, Haft DH, DiCuccio M, Chetvernin V, Badretdin A, Coulouris G, Chitsaz F, Derbyshire MK, Durkin AS, Gonzales NR, Gwadz M, Lanczycki CJ, Song JS, Thanki N, Wang J, Yamashita RA, Yang M, Zheng C, Marchler-Bauer A, Thibaud-Nissen F (2021). RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. Nucleic Acids Res. Jan 8;49 (D1):D1020-D1028. doi: 10.1093/nar/gkaa1105. PMID: 33270901; PMCID: PMC7779008.
"
input_fasta="$input_name"_"$database"_non_match.faa
if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_nfam_out "$DIR"/db/nfam/NF.hmm "$input_fasta"
else
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_nfam_out "$DIR"/db/nfam/NF.hmm "$input_fasta" > log_hmm_nfam
fi
check_success
hmmscan_out=$(<hmmscan_nfam_out)
database=nfam
map_file=NFAM_EC3.tsv
cut_column=2
hmm_func_annot
match=$(echo "$hmmscan_out" | list_hmm_hits)
report
echo "$m_count match(es) with NFAM" >> arche_report


# TIGRFAM database.
#-----------------

echo "

${bold}Queriying TIGRFAMs with HMMER3${reset} — Li W, O'Neill KR, Haft DH, DiCuccio M, Chetvernin V, Badretdin A, Coulouris G, Chitsaz F, Derbyshire MK, Durkin AS, Gonzales NR, Gwadz M, Lanczycki CJ, Song JS, Thanki N, Wang J, Yamashita RA, Yang M, Zheng C, Marchler-Bauer A, Thibaud-Nissen F (2021). RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. Nucleic Acids Res. Jan 8;49 (D1):D1020-D1028. doi: 10.1093/nar/gkaa1105. PMID: 33270901; PMCID: PMC7779008.
"
input_fasta="$input_name"_"$database"_non_match.faa
if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_tigrfam_out "$DIR"/db/tigrfam/TIGRFAM.hmm "$input_fasta"
else
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_tigrfam_out "$DIR"/db/tigrfam/TIGRFAM.hmm "$input_fasta" > log_hmm_tigrfam
fi
check_success
hmmscan_out=$(<hmmscan_tigrfam_out)
database=tigrfam
map_file=TIGRFAM_EC3.tsv
cut_column=2
hmm_func_annot
match=$(echo "$hmmscan_out" | list_hmm_hits)
report
echo "$m_count match(es) with TIGRFAM" >> arche_report


# Pfam database.
#-------------

echo "

${bold}Queriying PFAM with HMMER3${reset} — J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G.A. Salazar, E.L.L. Sonnhammer, S.C.E. Tosatto, L. Paladin, S. Raj, L.J. Richardson, R.D. Finn, A. Bateman. (2020) Pfam: The protein families database in 2021. Nucleic Acids Research. doi: 10.1093/nar/gkaa913
"
input_fasta="$input_name"_"$database"_non_match.faa
if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_pfam_out "$DIR"/db/pfam/pfam-A.hmm "$input_fasta"
else
	hmmscan --noali --cut_ga --cpu "$cpus" --domtblout hmmscan_pfam_out "$DIR"/db/pfam/pfam-A.hmm "$input_fasta" > log_hmm_pfam
fi
check_success
hmmscan_out=$(<hmmscan_pfam_out)
database=pfam
map_file=PFAM_EC7.tsv
cut_column=2
hmm_func_annot
match=$(echo "$hmmscan_out" | list_hmm_hits)
report
echo "$m_count match(es) with PFAM" >> arche_report
final_non_match=$(cat "$input_name"_"$database"_non_match.faa | grep ">" | wc -l)
echo "

${bold}${yellow}$final_non_match features didn't matched to any database and will be labelled as 'Hypothetical protein'.${reset}
"
echo "$final_non_match hypotetical protein(s)" >> arche_report

# Summary of tables.
#----------------

echo "${bold}Summarizing tables ... ${reset}
"
grep ">" "$input_name"_"$database"_non_match.faa | cut -f1 -d " " | sed 's/>//g ; s/$/| - | Hypothetical protein |/g' > non_match_TABLE
cat *_TABLE > summ_TABLE
locus_list=$(grep ">" "$input_name"_struc_annot.faa | sed 's/>//g' | cut -f1 -d " ")

if [[ $annotator = genemark ]]; then
	grep ">" "$input_name"_struc_annot.faa | cut -f1-4,6 -d " " | sed 's/>//g ; s/gene_type=/|/g ; s/ /|/ ; s/$/|/g' > pre_ORDER
	echo "$locus_list" | while read -r line ; do grep -w -m 1 "$line" summ_TABLE >> post_ORDER || printf "\n" >> post_ORDER ; done
	paste pre_ORDER post_ORDER | cut -f2-3,5-12 -d "|" > prot_FINAL
else	
	echo "$locus_list" | while read -r line ; do grep -w -m 1 "$line" summ_TABLE >> prot_FINAL || printf "\n" >> prot_FINAL ; done
fi

if [[ $annotator = genemark ]]; then
	grep ">" tRNA.fna | cut -f2,4-5 -d " " | sed 's/ /|tRNA / ; s/:/ /g ; s/-/ /g; s/$/| - |Transfer RNA|/g' > tRNA_FINAL 
	paste <(grep ">" rRNA.fna | cut -f3,4 -d ":" | sed 's/(.*//g; s/$/|/g ; s/:/ /g ; s/-/ /g') <(grep ">" rRNA.fna | cut -f1 -d ":" | sed 's/>//g ; s/_/ /g; s/$/| - |Ribosomal RNA|/g') > rRNA_FINAL
	output=$(cat *_FINAL | sed 's/ /:/; s/ /-/' | sed -E 's/\t//g' | sort -V | nl )
else
	output=$(cat *_FINAL | sed -E 's/\t//g' | sort -V | nl )
fi



echo "$output" > "$input_name"_omic_table.tbl
echo "$output" | sed 's/|/\t/g' > "$input_name"_omic_table.tsv
rm *_ORDER *_TABLE *_FINAL "${!#}" &>/dev/null

echo "Work completed.

Output files had been saved in "$directory_name"_"$random_append"

Final tables with annotations include ${bold}${green}"$input_name"_omic_table.tbl${reset} and ${bold}${green}"$input_name"_omic_table.tsv ${reset}
"
end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "
Total duration of the process: $elapsed seconds.
"
echo "Total duration of the process: $elapsed seconds." >> arche_report

cd "$current_dir"
