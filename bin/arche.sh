#!/bin/bash

start_time="$(date -u +%s)"

#Default parameters
#---------------

# Main directory
DIR=/home/cisme/Documentos/arche_1.0.1
# Procedure
mode=standard
kingdom=bact
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


#-------------------------------------------------------------------------------------------------------------------------------------------------
#Functions.
#-------------------------------------------------------------------------------------------------------------------------------------------------
usage() {
echo "
arche 1.0
---------
Written by Daniel Alonso <gundizalvus16@hotmail.com>
Last modified November 10, 2022
Description:
        Optimized annotator for prokariotic meta(genomes).

Usage:
	arche.sh --install [stand|compl]
        arche.sh --help
	arche.sh [options] <contigs.fasta>

Options:
	-h, --help           This help.
	-i, --install	     Set up the executable location, and install databases.
	-n, --name-files     Name of the files to be created in the output directory, in-
			     cluding the directory itself (default 'arche').
	-o, --output	     Provide the full path to the directory where the output di-
			     rectory will be created. E.g. /home/user/ (default current).
	-k, --kingdom        Source of the contigs. Use 'arch' for archaeal genomes or
                             'meta' for metagenomes (default is for bacterial genomes).
	-m, --mode           Gives priority to Orthology (KO, eggNOG) or Enzyme Comission
	                     designed databases during the annotation. Use 'kegg' for KO-->
                             eggNOG-->E.C., 'eggnog' for eggNOG-->KO-->E.C., or 'ec' for
                             E.C.-->KO-->eggNOG (default will use a shorter swiss-prot KO·
                             ·eggNOG·E.C. designed database with no priority).
	-a, --alignment      Select the algorithm to use during the protein alignment step:
                             'diamond' (accelerated blastp) or 'ssearch' (Smith-Waterman)
                             (default 'blastp').
	-t, --threads        Number of threads to use (default '1').
	-r, --memory         Amount of RAM to use in GB (default '2').
	-e, --evalue         Similarity e-value cut-off (default '1e-08').
	-q, --query-cov      Minimum coverage on query protein (default '70').
	-b, --bypass         Bypass the RNA gene prediction (default 'off').
	-v, --verbose	     Verbose mode (default 'off')
	"
}

install_dir () {
main_dir=$(realpath arche.sh | sed 's!/bin/arche.sh!!')
	sed -i -z "s|DIR=[^\n]*|DIR=$main_dir|" arche.sh
}

install_db () {
echo "
${bold}${cyan}Unpacking databases and map files...${reset}
"
7z e ../db/nfam/NF.hmm.7z -o../db/nfam/
rm ../db/nfam/NF.hmm.7z
7z e ../db/ssearch/EC_db.fa.7z -o../db/ssearch/
7z e ../db/ssearch/EC_db_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/EC_db_MET.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG-EC_db.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG-EC_db_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG-EC_db_MET.fa.7z -o../db/ssearch/
7z e ../db/ssearch/eggNOG-EC_db.fa.7z -o../db/ssearch/
7z e ../db/ssearch/eggNOG-EC_db_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/eggNOG-EC_db_MET.fa.7z -o../db/ssearch/
7z e ../db/ssearch/sprot.fa.7z -o../db/ssearch/
7z e ../db/ssearch/sprot_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/sprot_MET.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG_db.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG_db_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG_db_MET.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG-eggNOG_db.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG-eggNOG_db_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/KEGG-eggNOG_db_MET.fa.7z -o../db/ssearch/
7z e ../db/ssearch/eggNOG_db.fa.7z -o../db/ssearch/
7z e ../db/ssearch/eggNOG_db_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/eggNOG_db_MET.fa.7z -o../db/ssearch/
7z e ../db/ssearch/EC-KEGG-eggNOG_db.fa.7z -o../db/ssearch/
7z e ../db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa.7z -o../db/ssearch/
7z e ../db/ssearch/EC-KEGG-eggNOG_db_MET.fa.7z -o../db/ssearch/
rm ../db/ssearch/*.7z
7z e ../db/kofam/KO_prokaryote.hmm.7z -o../db/kofam/
rm ../db/kofam/KO_prokaryote.hmm.7z
7z e ../db/tigrfam/TIGRFAM.hmm.7z -o../db/tigrfam/
rm ../db/tigrfam/TIGRFAM.hmm.7z
7z e ../db/pfam/pfam-A.hmm.7z -o../db/pfam/
rm ../db/pfam/pfam-A.hmm.7z
7z e ../mpp/NFAM_EC3.tsv.7z -o../mpp/
7z e ../mpp/uniprot_bac.7z -o../mpp/
7z e ../mpp/uniprot_arch.7z -o../mpp/
7z e ../mpp/uniprot_met.7z -o../mpp/
7z e ../mpp/KO_match_table6.tsv.7z -o../mpp/
7z e ../mpp/TIGRFAM_EC3.tsv.7z -o../mpp/
7z e ../mpp/PFAM_EC7.tsv.7z -o../mpp/
rm ../mpp/*.7z

echo "
${bold}${cyan}Setting up databases...${reset}
"
diamond makedb --in ../db/ssearch/sprot.fa -d ../db/diamond/sprot
diamond makedb --in ../db/ssearch/sprot_ARCH.fa -d ../db/diamond/sprot_ARCH
diamond makedb --in ../db/ssearch/sprot_MET.fa -d ../db/diamond/sprot_MET
diamond makedb --in ../db/ssearch/EC_db.fa -d ../db/diamond/EC_db
diamond makedb --in ../db/ssearch/EC_db_ARCH.fa -d ../db/diamond/EC_db_ARCH
diamond makedb --in ../db/ssearch/EC_db_MET.fa -d ../db/diamond/EC_db_MET
diamond makedb --in ../db/ssearch/KEGG_db.fa -d ../db/diamond/KEGG_db
diamond makedb --in ../db/ssearch/KEGG_db_ARCH.fa -d ../db/diamond/KEGG_db_ARCH
diamond makedb --in ../db/ssearch/KEGG_db_MET.fa -d ../db/diamond/KEGG_db_MET
diamond makedb --in ../db/ssearch/eggNOG_db.fa -d ../db/diamond/eggNOG_db
diamond makedb --in ../db/ssearch/eggNOG_db_ARCH.fa -d ../db/diamond/eggNOG_db_ARCH
diamond makedb --in ../db/ssearch/eggNOG_db_MET.fa -d ../db/diamond/eggNOG_db_MET
diamond makedb --in ../db/ssearch/KEGG-eggNOG_db.fa -d ../db/diamond/KEGG-eggNOG_db
diamond makedb --in ../db/ssearch/KEGG-eggNOG_db_ARCH.fa -d ../db/diamond/KEGG-eggNOG_db_ARCH
diamond makedb --in ../db/ssearch/KEGG-eggNOG_db_MET.fa -d ../db/diamond/KEGG-eggNOG_db_MET
diamond makedb --in ../db/ssearch/KEGG-EC_db.fa -d ../db/diamond/KEGG-EC_db
diamond makedb --in ../db/ssearch/KEGG-EC_db_ARCH.fa -d ../db/diamond/KEGG-EC_db_ARCH
diamond makedb --in ../db/ssearch/KEGG-EC_db_MET.fa -d ../db/diamond/KEGG-EC_db_MET
diamond makedb --in ../db/ssearch/eggNOG-EC_db.fa -d ../db/diamond/eggNOG-EC_db
diamond makedb --in ../db/ssearch/eggNOG-EC_db_ARCH.fa -d ../db/diamond/eggNOG-EC_db_ARCH
diamond makedb --in ../db/ssearch/eggNOG-EC_db_MET.fa -d ../db/diamond/eggNOG-EC_db_MET
diamond makedb --in ../db/ssearch/EC-KEGG-eggNOG_db.fa -d ../db/diamond/EC-KEGG-eggNOG_db
diamond makedb --in ../db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa -d ../db/diamond/EC-KEGG-eggNOG_db_ARCH
diamond makedb --in ../db/ssearch/EC-KEGG-eggNOG_db_MET.fa -d ../db/diamond/EC-KEGG-eggNOG_db_MET

makeblastdb -in ../db/ssearch/sprot.fa -dbtype prot -parse_seqids -out ../db/blastp/sprot
makeblastdb -in ../db/ssearch/sprot_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/sprot_ARCH
makeblastdb -in ../db/ssearch/sprot_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/sprot_MET
makeblastdb -in ../db/ssearch/EC_db.fa -dbtype prot -parse_seqids -out ../db/blastp/EC_db
makeblastdb -in ../db/ssearch/EC_db_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/EC_db_ARCH
makeblastdb -in ../db/ssearch/EC_db_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/EC_db_MET
makeblastdb -in ../db/ssearch/KEGG_db.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG_db
makeblastdb -in ../db/ssearch/KEGG_db_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG_db_ARCH
makeblastdb -in ../db/ssearch/KEGG_db_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG_db_MET
makeblastdb -in ../db/ssearch/eggNOG_db.fa -dbtype prot -parse_seqids -out ../db/blastp/eggNOG_db
makeblastdb -in ../db/ssearch/eggNOG_db_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/eggNOG_db_ARCH
makeblastdb -in ../db/ssearch/eggNOG_db_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/eggNOG_db_MET
makeblastdb -in ../db/ssearch/KEGG-eggNOG_db.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG-eggNOG_db
makeblastdb -in ../db/ssearch/KEGG-eggNOG_db_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG-eggNOG_db_ARCH
makeblastdb -in ../db/ssearch/KEGG-eggNOG_db_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG-eggNOG_db_MET
makeblastdb -in ../db/ssearch/KEGG-EC_db.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG-EC_db
makeblastdb -in ../db/ssearch/KEGG-EC_db_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG-EC_db_ARCH
makeblastdb -in ../db/ssearch/KEGG-EC_db_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/KEGG-EC_db_MET
makeblastdb -in ../db/ssearch/eggNOG-EC_db.fa -dbtype prot -parse_seqids -out ../db/blastp/eggNOG-EC_db
makeblastdb -in ../db/ssearch/eggNOG-EC_db_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/eggNOG-EC_db_ARCH
makeblastdb -in ../db/ssearch/eggNOG-EC_db_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/eggNOG-EC_db_MET
makeblastdb -in ../db/ssearch/EC-KEGG-eggNOG_db.fa -dbtype prot -parse_seqids -out ../db/blastp/EC-KEGG-eggNOG_db
makeblastdb -in ../db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa -dbtype prot -parse_seqids -out ../db/blastp/EC-KEGG-eggNOG_db_ARCH
makeblastdb -in ../db/ssearch/EC-KEGG-eggNOG_db_MET.fa -dbtype prot -parse_seqids -out ../db/blastp/EC-KEGG-eggNOG_db_MET

map_db ../db/ssearch/sprot.fa
map_db ../db/ssearch/sprot_ARCH.fa
map_db ../db/ssearch/sprot_MET.fa
map_db ../db/ssearch/EC_db.fa
map_db ../db/ssearch/EC_db_ARCH.fa
map_db ../db/ssearch/EC_db_MET.fa
map_db ../db/ssearch/KEGG_db.fa
map_db ../db/ssearch/KEGG_db_ARCH.fa
map_db ../db/ssearch/KEGG_db_MET.fa
map_db ../db/ssearch/eggNOG_db.fa
map_db ../db/ssearch/eggNOG_db_ARCH.fa
map_db ../db/ssearch/eggNOG_db_MET.fa
map_db ../db/ssearch/KEGG-eggNOG_db.fa
map_db ../db/ssearch/KEGG-eggNOG_db_ARCH.fa
map_db ../db/ssearch/KEGG-eggNOG_db_MET.fa
map_db ../db/ssearch/KEGG-EC_db.fa
map_db ../db/ssearch/KEGG-EC_db_ARCH.fa
map_db ../db/ssearch/KEGG-EC_db_MET.fa
map_db ../db/ssearch/eggNOG-EC_db.fa
map_db ../db/ssearch/eggNOG-EC_db_ARCH.fa
map_db ../db/ssearch/eggNOG-EC_db_MET.fa
map_db ../db/ssearch/EC-KEGG-eggNOG_db.fa
map_db ../db/ssearch/EC-KEGG-eggNOG_db_ARCH.fa
map_db ../db/ssearch/EC-KEGG-eggNOG_db_MET.fa

hmmpress ../db/kofam/KO_prokaryote.hmm
hmmpress ../db/pfam/pfam-A.hmm
hmmpress ../db/tigrfam/TIGRFAM.hmm
hmmpress ../db/nfam/NF.hmm
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
		[ $? -eq 0 ] || { 
		echo "ERROR: blastp. Check options!" ; exit 1
		}
		match=$(awk -v cov2=$cov '{if ($3 > cov2) print $1}' "$database"_out | grep -v "#")
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
		[ $? -eq 0 ] || { 
		echo "ERROR: blastp. Check options!" ; exit 1 
		}
		match=$(cat "$database"_out | grep -v "*" | cut -f1)
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
		[ $? -eq 0 ] || exit 1
		match=$(grep -v "#" "$database"_out | awk -v cov2=$cov '{if ((($6/$2)*100 > cov2)) print $1"\t"$3}' | cut -f1)	
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
		echo -n "${bold}${cyan}0 match(es) reported${reset}"
		cp "$input_fasta" "$input_name"_"$database"_non_match.faa
	fi
}

cite_blastp () {
	echo "
---------  ${bold}Querying Uniprot with BLASTp+ (heuristic algorithm)${reset} ----------
	
The UniProt Consortium. UniProt: the universal protein knowledgebase in 2021 (2021). Nuc-
leic Acids Res. 49:D1
	
Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. 
(2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421.
"
}

cite_diamond () { 
	echo "
	
------- ${bold}Queriying Uniprot with DIAMOND (fast heuristic algorithm)${reset}  -------
	
The UniProt Consortium. UniProt: the universal protein knowledgebase in 2021 (2021). Nucleic 
Acids Res. 49:D1

Buchfink B, Xie C, Huson DH, 'Fast and sensitive protein alignment using DIAMOND', Nature Me-
thods 12, 59-60 (2015). doi:10.1038/nmeth.3176
"
}

cite_ssearch () {
echo "
------  ${bold}Queriying Uniprot with SSEARCH (Smith-Waterman algorithm)${reset}  ------
	
The UniProt Consortium. UniProt: the universal protein knowledgebase in 2021 (2021). Nucleic 
Acids Res. 49:D1
	
W. R. Pearson. Effective protein sequence comparison (1996). Methods Enzymol., 266:227–258
"
}

kofam_search () {
	echo "
	
-------------  ${bold}Querying KofamKOALA (prokaryotic) with HMMER3${reset}  -----------

S. R. Eddy. (2011)  Accelerated profile HMM searches. PLoS Comp. Biol., 7:e1002195.

Takuya Aramaki, Romain Blanc-Mathieu, Hisashi Endo, Koichi Ohkubo, Minoru Kanehisa, Susumu 
Goto, Hiroyuki Ogata. (2020) KofamKOALA: KEGG Ortholog assignment based on profile HMM and 
adaptive score threshold. Bioinformatics. 36:2251–2252
"
	input_fasta="$input_name"_"$database"_non_match.faa
	if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_kofam_out "$DIR"/db/kofam/KO_prokaryote.hmm "$input_fasta"
	else
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_kofam_out "$DIR"/db/kofam/KO_prokaryote.hmm "$input_fasta" > log_hmm_kofam	
	fi
	[ $? -eq 0 ] || exit 1
	hmmscan_out=$(<hmmscan_kofam_out)
	database=kofam
	map_file=KO_match_table6.tsv
	cut_column=1
	hmm_func_annot
	match=$(echo "$hmmscan_out" | list_hmm_hits)
	report ; echo "${bold}${cyan} with KEGG only associated number(s) (from HMMs).${reset}"
}

hmm_func_annot () {
	echo "$hmmscan_out" | sed '/^#/d' | awk '!x[$4]++' | awk -v cov3=$cov -v column=$cut_column '{if ((($17-$16)/$6)*100 > cov3) print $column"\t"$4}' | sed -E 's/\t/|/g' |  sed 's/$/|/g' | sort  > pre_TABLE
	map_LIST=$(cut -f1 -d "|" pre_TABLE | cut -f1 -d ".")
	echo "$map_LIST" | while read -r line ; do grep -w -m 1 "$line" "$DIR"/mpp/"$map_file" >> post_TABLE || printf "\n" >> post_TABLE ; done
	paste pre_TABLE post_TABLE  | sed -E 's/\t//g' | cut -d "|" -f2-10 > "$database"_TABLE
	rm pre_TABLE post_TABLE
}

list_hmm_hits () {
	sed '/^#/d' | awk '!x[$4]++' | awk -v cov3=$cov '{if ((($17-$16)/$6)*100 > cov3) print $1"\t"$4}' | cut -f2
}



#-------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------


#Options
#------

options=$(getopt -l "install,kingdom:,bypass,mode;,alignment:,name-files:,output:,threads:,evalue:,memory:,query-cov:,help,verbose" -o "ik:bm:a:n:o:t:e:r:q:hv" -- "$@")
[ $? -eq 0 ] || {
	echo "Incorrect options provided."
    exit 1
}
eval set -- "$options"
while true ; do
	case $1 in
	-i|--install)
		shift; install_dir; install_db
		exit 0
		;;
	-k|--kingdom)
		shift
		[[ ! "$1" =~ ^arch$|^meta$ ]] && {
            echo "Incorrect option provided for '-k'/'--kingdom'. Choose 'arch' for archaea or 'meta' for metagenome. Remove the option for bacterial genomes."
            exit 1
        }
		 kingdom=$1
		;;
	-b|--bypass)
		shift; bypass="true"
		;;
	-m|--mode)
		shift
		[[ ! "$1" =~ ^kegg$|^eggnog$|^ec$ ]] && {
            echo "Incorrect option provided for '-m'/'--mode'. Choose 'kegg', 'eggnog' or 'ec'."
            exit 1
        }
		 mode=$1
		;;
	-a|--alignment)
		shift
		[[ ! "$1" =~ ^diamond$|^ssearch$ ]] && {
            echo "Incorrect option provided for '-a'/'--engine'. Choose 'diamond' or 'ssearch'."
            exit 1
        }
		 db_search_alig=$1
		;;
	-n|--name-files)
		shift; input_name=$1
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
		shift; cpus=$1
		;;
	-e|--evalue)
		shift; evalue=$1
		;;
	-r|--memory)
		shift; memory=$1
		;;
	-q|--query-cov)
		shift; cov=$1
		;;
	-h|--help)
		shift; usage
		exit 0
		;;
	-v|--verbose)
		shift; 
		quiet_on=
		quiet_on_prodigal=
		verbose_on=-v
		;;
	--)
		shift
	    break;;
	esac
	shift
done



# Install-configure
#--------------
	


# In case of archaeal genome, or metagenome
#------------------------------------

if [[ "$kingdom" = arch ]]; then	
	 barrnap_mode=arc
	 tRNA_mode=-A
	 gms2_mode=archaea
	uniprot_db=arch
	uni_db_kingdom=_ARCH
elif [[ "$kingdom" = meta ]]; then	
	 annotator=metaprodigal
	 db_search_alig=diamond
	uniprot_db=met
	uni_db_kingdom=_MET
fi


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
${bold}Creating output directory "$directory_name"_"$random_append"${reset}
"
	mkdir  "$dir_path""$directory_name"_"$random_append"
	cp "${!#}" "$dir_path""$directory_name"_"$random_append"/
	cd "$dir_path""$directory_name"_"$random_append"
fi
	

# Check contig name for further "sort -V"
#----------------------------------
if [[ "$kingdom" != meta ]]; then
	contig_file="${!#}"
	contig_name_start=$(grep ">" "$contig_file" | head -n 1 | sed 's/>//' | cut -c1);
	if [[ "$contig_name_start" != [A-Za-z] ]]; then
		echo "
${bold}${yellow}Contig headers start with non-letter characters... creating a new contig file ("$contig_file"_2). ${reset}"
		sed 's/>/>c_/g' "$contig_file" > "$contig_file"_2
		# [ $? -eq 0 ] || exit 1
		contig_file="$contig_file"_2
	fi
fi




if [[ "$bypass" != "true" ]]; then

# Ribosomal RNA gene prediction.
#--------------------------

	[ $annotator = genemark ] && {	
		echo "

----------------  ${bold}Predicting ribosomal RNA with Barrnap${reset}  ----------------

Seemann T. barrnap 0.9: rapid ribosomal RNA prediction https://github.com/tseemann/barrnap
"
		barrnap $quiet_on --kingdom  "$barrnap_mode" --outseq rRNA.fna --threads "$cpus" "$contig_file" > rRNA.tsv
		[ $? -eq 0 ] || { 
		echo "ERROR: barrnap." ; exit 1 
		}
		barrnap_exit=$(grep ">" rRNA.fna | wc -l)
		echo "
${bold}${cyan}$barrnap_exit rRNAs gene(s) reported.${reset}"
}	



# tRNA gene prediction.
#------------------

	[ $annotator = genemark ] && {
		echo "
---------------  ${bold}Predicting transfer RNA with tRNAscan-SE${reset}  --------------

Chan, P.P. and Lowe, T. M. (2019) tRNAscan-SE: Searching for tRNA Genes in Genomic Seque-
nces. Methods Mol Biol. 1962:1-14.
"
		tRNAscan-SE $tRNA_mode $quiet_on -o tRNA.tsv -a tRNA.fna "$contig_file"
		[ $? -eq 0 ] || exit 1
		tRNAscan_exit=$(grep ">" tRNA.fna | wc -l)
		echo "
${bold}${cyan}$tRNAscan_exit tRNAs gene(s) reported.${reset}"
}

fi




# Structural annotation.
#------------------

if [[ $annotator = genemark ]]; then 	
	echo "
---------------  ${bold}Structural annotation with GeneMarkS-2${reset}  ----------------

Lomsadze A, Gemayel K, Tang S, Borodovsky M (2018). Modeling leaderless transcription and 
atypical genes results in more accurate gene prediction in prokaryotes. Genome Res, 29(7), 
pp 1079-1089
"	
	gms2.pl --gcode 11 --seq  "$contig_file" --fnn "$input_name"_struc_annot.fna --faa "$input_name"_struc_annot.faa  --genome-type "$gms2_mode" --output "$input_name"_GeneMark.lst
	[ $? -eq 0 ] || exit 1
	features_tot=$(cat "$input_name"_struc_annot.fna | grep ">" | wc -l)
	atypical_tot=$(cat "$input_name"_struc_annot.fna | grep "atypical" | wc -l)
	let native=$features_tot-$atypical_tot
	echo "${bold}${cyan}$features_tot feature(s) reported; $native were native and $atypical_tot atypical.${reset}"
else
	echo "

---------------  ${bold}Structural annotation with MetaProdigal ${reset}  ----------------

Hyatt,D., LoCascio,P.F., Hauser,L.J. and Uberbacher,E.C. (2012) Gene and translation initiation
site prediction in metagenomic sequences. Bioinformatics, 28, 2223–2230
" 
	prodigal -g 11 -a "$input_name"_struc_annot.faa -d "$input_name"_struc_annot.fna -f gff -o "$input_name"_struc_annot.gff $quiet_on_prodigal -p meta -i "${!#}" > prodigal.log
	[ $? -eq 0 ] || exit 1
	features_tot=$(cat "$input_name"_struc_annot.fna | grep ">" | wc -l)
	echo "
${bold}${cyan}$features_tot feature(s) reported.${reset}"
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
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic1
	db_heurist=db01 #KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with KEGG + eggNOG associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic2
	db_heurist=db02 #KEGG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with KEGG + E.C. associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic3
	db_heurist=db03 #KEGG.tsv
	heurist_search ; echo "${bold}${cyan} with KEGG only associated number(s) (from alignments).${reset}"
	
	# Kofam database.
	kofam_search ; echo ""
	
	# ... continuation.
	cite_"$db_search_alig"
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic4
	db_heurist=db04 #eggNOG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + E.C. associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic5
	db_heurist=db05
	heurist_search ; echo "${bold}${cyan} with eggNOG only associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic6
	db_heurist=db06 #EC.tsv
	heurist_search ; echo "${bold}${cyan} with E.C. only associated number(s).${reset}"


	;;
eggnog)
	# Prioritize target sequences with eggNOG associated IDs.
	cite_"$db_search_alig"
	input_fasta="$input_name"_struc_annot.faa
	database=heuristic0
	db_heurist=db00 #EC-KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + KEGG + E.C. associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic1
	db_heurist=db01 #KEGG-eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + KEGG associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic2
	db_heurist=db04 #eggNOG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG + E.C. associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic3
	db_heurist=db05 #eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG only associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic4
	db_heurist=db02 #KEGG-EC
	heurist_search ; echo "${bold}${cyan} with KEGG + E.C. associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic5
	db_heurist=db03 #KEGG
	heurist_search ; echo "${bold}${cyan} with KEGG only associated number(s) (from alignments).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic6
	db_heurist=db06 #EC
	heurist_search ; echo "${bold}${cyan} with E.C. only associated number(s).${reset}"

	
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
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic1
	db_heurist=db02 #KEGG-EC
	heurist_search ; echo "${bold}${cyan} with E.C. + KEGG associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic2
	db_heurist=db04 #eggNOG-EC.tsv
	heurist_search ; echo "${bold}${cyan} with E.C. + eggNOG associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic3
	db_heurist=db06 #EC.tsv
	heurist_search ; echo "${bold}${cyan} with E.C. only associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic4
	db_heurist=db01 #KEGG-eggNOG
	heurist_search ; echo "${bold}${cyan} with KEGG + eggNOG associated number(s).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic5
	db_heurist=db03 #KEGG
	heurist_search ; echo "${bold}${cyan} with KEGG only associated number(s) (from alignments).${reset}"
	
	input_fasta="$input_name"_"$database"_non_match.faa
	database=heuristic6
	db_heurist=db05 #eggNOG.tsv
	heurist_search ; echo "${bold}${cyan} with eggNOG only associated number(s).${reset}"
	
	# Kofam database.
	kofam_search
	;;
	
esac



# NFAM database.
#--------------

echo "

----------- ${bold}Querying NCBI Protein Family Models with HMMER3${reset}  -----------

Li W, O'Neill KR, Haft DH, DiCuccio M, Chetvernin V, Badretdin A, Coulouris G, Chitsaz F,
Derbyshire MK, Durkin AS, Gonzales NR, Gwadz M, Lanczycki CJ, Song JS, Thanki N, Wa-
ng J, Yamashita RA, Yang M, Zheng C, Marchler-Bauer A, Thibaud-Nissen F (2021). Ref-
Seq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family mo-
del curation. Nucleic Acids Res. Jan 8;49 (D1):D1020-D1028. doi: 10.1093/nar/gkaa1105. 
PMID: 33270901; PMCID: PMC7779008.
"
input_fasta="$input_name"_"$database"_non_match.faa
if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_nfam_out "$DIR"/db/nfam/NF.hmm "$input_fasta"
else
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_nfam_out "$DIR"/db/nfam/NF.hmm "$input_fasta" > log_hmm_nfam
fi
[ $? -eq 0 ] || exit 1
hmmscan_out=$(<hmmscan_nfam_out)
database=nfam
map_file=NFAM_EC3.tsv
cut_column=2
hmm_func_annot
match=$(echo "$hmmscan_out" | list_hmm_hits)
report



# TIGRFAM database.
#-----------------

echo "

------------------- ${bold}Queriying TIGRFAMs with HMMER3${reset}  --------------------

Li W, O'Neill KR, Haft DH, DiCuccio M, Chetvernin V, Badretdin A, Coulouris G, Chitsaz F,
Derbyshire MK, Durkin AS, Gonzales NR, Gwadz M, Lanczycki CJ, Song JS, Thanki N, Wa-
ng J, Yamashita RA, Yang M, Zheng C, Marchler-Bauer A, Thibaud-Nissen F (2021). Ref-
Seq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family mo-
del curation. Nucleic Acids Res. Jan 8;49 (D1):D1020-D1028. doi: 10.1093/nar/gkaa1105. 
PMID: 33270901; PMCID: PMC7779008.
"
input_fasta="$input_name"_"$database"_non_match.faa
if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_tigrfam_out "$DIR"/db/tigrfam/TIGRFAM.hmm "$input_fasta"
else
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_tigrfam_out "$DIR"/db/tigrfam/TIGRFAM.hmm "$input_fasta" > log_hmm_tigrfam
fi
[ $? -eq 0 ] || exit 1
hmmscan_out=$(<hmmscan_tigrfam_out)
database=tigrfam
map_file=TIGRFAM_EC3.tsv
cut_column=2
hmm_func_annot
match=$(echo "$hmmscan_out" | list_hmm_hits)
report



# Pfam database.
#-------------

echo "

----------------------  ${bold}Queriying PFAM with HMMER3${reset}  ---------------------

J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G.A. Salazar, E.L.L. Sonnhammer, S.C.
E. Tosatto, L. Paladin, S. Raj, L.J. Richardson, R.D. Finn, A. Bateman. (2020) Pfam: The 
protein families database in 2021. Nucleic Acids Research. doi: 10.1093/nar/gkaa913
"
input_fasta="$input_name"_"$database"_non_match.faa
if [[ -n "$verbose_on" ]]; then
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_pfam_out "$DIR"/db/pfam/pfam-A.hmm "$input_fasta"
else
	hmmscan --noali --incE "$evalue" -E "$evalue" --cpu "$cpus" --domtblout hmmscan_pfam_out "$DIR"/db/pfam/pfam-A.hmm "$input_fasta" > log_hmm_pfam
fi
[ $? -eq 0 ] || exit 1
hmmscan_out=$(<hmmscan_pfam_out)
database=pfam
map_file=PFAM_EC7.tsv
cut_column=2
hmm_func_annot
match=$(echo "$hmmscan_out" | list_hmm_hits)
report
final_non_match=$(cat "$input_name"_"$database"_non_match.faa | grep ">" | wc -l)
echo "

${bold}${yellow}$final_non_match features didn't matched to any database and will be labelled as 'Hypothetical protein'.${reset}
"


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
# rm *_ORDER *_TABLE *_FINAL "${!#}" &>/dev/null

echo "Work completed.

Output files had been saved in "$directory_name"_"$random_append"

Final tables with annotations include ${bold}${green}"$input_name"_omic_table.tbl${reset} and ${bold}${green}"$input_name"_omic_table.tsv ${reset}
"
cd "$current_dir"

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "
Total duration of the process: $elapsed seconds.
"

# exec bash
# rm -v !("arche.sh"|"pro_arche.sh"|"Nesterenkonia_contigs.fna"|"sub_arche.sh"|"arche_config"|"compost_contigs.fna"|"KF577590.1.fasta"|"Nesterenkonia_contigs1+2.fna"|"gut_meta2.fna"|"gut_meta.fna") ; exec bash
