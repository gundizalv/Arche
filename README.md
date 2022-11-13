[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI:10.1093/bioinformatics/btu153](https://zenodo.org/badge/DOI/10.1093/bioinformatics/btu153.svg)](https://doi.org/10.1093/bioinformatics/btu153)
![Don't judge me](https://img.shields.io/badge/Language-Perl_5-steelblue.svg)

# ARCHE: a functional-optimized annotator for microbial meta(genomes)

## Installing dependencies

Before you download ARCHE (13Gb), make sure the following software are working properly on your computer:

p7zip  
bedtools >= 2.27.0  
barrnap  
Prodigal  
hmmer  
blastp  
DIAMOND 

FASTA36 ---> https://github.com/wrpearson/fasta36.git  

faSomeRecords ---> The version which works with ARCHE can be found here http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/  

Infernal 1.1.4 ---> Install manually!  
&nbsp;    **tar xvfz infernal-1.1.4.tar.gz**  
&nbsp;    **cd infernal-1.1.4/**  
&nbsp;    **./configure**  
&nbsp;    **make**  
&nbsp;    **sudo make install**  
&nbsp;    **sudo mv infernal-1.1.4 /opt/**  

tRNAscan-SE ---> Install manually!  
&nbsp;    **Download it from http://lowelab.ucsc.edu/tRNAscan-SE/**  
&nbsp;    **tar xvfz trnascan-se-2.0.9.tar.gz**  
&nbsp;    **cd tRNAscan-SE-2.0**  
&nbsp;    **make**  
&nbsp;    **sudo make install**  
&nbsp;    **sudo mv tRNAscan-SE-2.0/ /opt/**  

GeneMarkS-2  
&nbsp;    **Download GeneMarkS-2 and key from http://exon.gatech.edu/GeneMark/license_download.cgi**  
&nbsp;    **gunzip gm_key_64.gz**  
&nbsp;    **tar xvfz gms2_linux_64.tar.gz**  
&nbsp;    **cp gm_key_64 ~/.gmhmmp2_key**  
&nbsp;    **sudo mv gms2_linux_64 /opt/**  
&nbsp;    **cp -s /opt/gms2_linux_64/{biogem,compp,gmhmmp2,gms2.pl} /usr/bin/**  

## Running ARCHE

### BlastP annotation of a bacterial genome, using 20 threads and 40 GB of memory:
```
arche.sh -n ecoli -t 20 -r 40 e_coli.fna
```

### SSEARCH annotation of an archaeal genome, using 1 thread and 2 GB of memory
```
arche.sh -n halorubrum -a ssearch -k achaea halorubrum_sp_DM2.fa
```

### DIAMOND annotation of metagenome
```
arche.sh -n seawater_meatgenome -k meta seawater_metagenome.fna
```


## Output Files

| File(s) | Description |
| --------- | ----------- |
| rRNA.tsv | GFF v3 file containing rRNA annotations. |
| rRNA.fna | FASTA file of all rRNA features. |
| tRNA.tsv | Table with tRNA details (coordinates, isotype, anticodon, scores, etc). |
| []_struc_annot.fna | FASTA file of all genomic features (nucleotide). |
| []_struc_annot.faa | FASTA file of translated coding genes (aminoacid). |
| heuristic[]_out | Output matches of the search instance(s) performed with BLASTp, DIAMOND or SSEARCH36. |
| heuristic[]_non_match.faa | FASTA file with the remaining non-matched sequences after the search instance(s) performed with BLASTp, DIAMOND or SSEARCH36. |
| hmmscan_[]_out | HMMER3 output table of the search instance(s) performed against a specific HMMDB. |
| [HMMDB]_non_match.faa | FASTA file with the remaining non-matched sequences after the search instance performed against a specific HMMDB. |
| []_omic_table.tbl | Feature table with fields separated by vertical bars. |
| []_omic_table.tsv | Feature table with tab-separated fields. |

## Command line options

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
    
	-b, --bypass         Use 'yes' to bypass the RNA gene prediction.
    
	-v, --verbose	     Verbose mode (default 'off')


# Licence

[GPL v3](https://raw.githubusercontent.com/tseemann/prokka/master/doc/LICENSE.Prokka)

## Author

* Daniel Alonso
