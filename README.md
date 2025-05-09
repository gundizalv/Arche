[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Don't judge me](https://img.shields.io/badge/Language-Bash-blue)
[![DOI:10.1101/2022.11.28.518280/bioRxiv](https://zenodo.org/badge/DOI/10.1101/2022.11.28.518280/bioRxiv.svg)](https://doi.org/10.1101/2022.11.28.518280)

# Arche: a flexible tool for annotation of microbial contigs

The growing amount of genomic data has prompted a need for less demanding and user friendly functional annotators. At the present, it’s hard to find a pipeline for the annotation of multiple functional data, such as both enzyme commission numbers (E.C.) and orthologous identifiers (KEGG and eggNOG), protein names, gene names, alternative names, and descriptions. Here, we provide a new solution which combines different algorithms (BLAST, DIAMOND, HMMER3) and databases (UniprotKB, KOfam, NCBIFAMs, TIGRFAMs, and PFAM), and also overcome data download challenges. Arche analysis pipeline can accommodate advanced tools in a unique order, creating several advantages regarding to other commonly used annotators.

## Installing GeneMarkS-2

Before you download Arche and its databases (13Gb), make sure GeneMarkS-2 (GMS2) is working properly on your computer. As GMS2 requires a licence (free), you must download it manually

&nbsp;    **Download GeneMarkS-2 and key from http://exon.gatech.edu/GeneMark/license_download.cgi**  

```
tar xvfz gms2_linux_[version].tar.gz
```
Move the dir to the desired place, and make the binary files accesible to your PATH.

Configure the key you've downloaded

```
gunzip gm_key.gz
cp gm_key ~/.gmhmmp2_key
``` 
or
```
cp gm_key ~/.gm_key
``` 

Test the software

```
gms2.pl --seq YOUR_GENOME
```

To install Arche, you will require the anaconda distribution. Download and install it from https://www.anaconda.com/download/success

```
conda create -n arche_env -c gundizalv92 -c bioconda -c conda-forge arche
conda activate arche_env
```


This command wil create a conda environment for arche future runs. It includes the installation of Arche from my user channel and other specific packages from bioconda and conda-forge channels.

## Installing databases

The already formatted databases and mapping files can be downloaded (13Gb) via command line:  

```
arche --download
```
Once you have the decompressed database directory (e.g. arche_0.0.0) do
```
arche --install
```


## Troubleshooting

In the case the instalation process or the running fails:

1. Check you are working within the conda environment you've created ("conda activate arche_env")
2. Check you have properly installed GeneMarkS-2

## Running Arche

### BlastP annotation of a bacterial genome, using 20 threads and 40 GB of memory:
```
arche -n ecoli -t 20 -r 40 e_coli.fna
```

### SSEARCH annotation of an archaeal genome, using 1 thread and 2 GB of memory
```
arche -n halorubrum -a ssearch -k achaea halorubrum_sp_DM2.fa
```

### DIAMOND annotation of a metagenome
```
arche -n seawater_meatgenome -k meta seawater_metagenome.fna
```
## Annotation of Escherichia coli K12

Here you can download a sample which includes the annotation of Escherichia coli K12 with several tools including Arche:

[https://docs.google.com/spreadsheets/d/17Nd_y7w2axfxsjFJYAvb_NI3AW9HjNx4/edit?usp=sharing&ouid=115908476093915484477&rtpof=true&sd=true](https://docs.google.com/spreadsheets/d/1ISnksbKhaqlUVYyr8r4OAOWdpsfri2F9/edit?usp=sharing&ouid=115908476093915484477&rtpof=true&sd=true)

## Output Table

The final table with all the annotations comes in two flavours:

[...]_omic_table.tbl which can be examined through the linux console with the comand 

&nbsp;    **column -ts "|" [example]_omic_table.tbl | less -S**  

 [...]_omic_table.tsv which can be opened using spreadsheet editors like Microsoft Excel, LibreOffice Calc, etc.


## Output Files

| File(s) | Description |
| --------- | ----------- |
| rRNA.tsv | GFF v3 file containing rRNA annotations. |
| rRNA.fna | FASTA file of all rRNA features. |
| tRNA.tsv | Table with tRNA details (coordinates, isotype, anticodon, scores, etc). |
| [...]_struc_annot.fna | FASTA file of all genomic features (nucleotide). |
| [...]_struc_annot.faa | FASTA file of translated coding genes (aminoacid). |
| heuristic[...]_out | Output matches of the search instance(s) performed with BLASTp, DIAMOND or SSEARCH36. |
| heuristic[...]_non_match.faa | FASTA file with the remaining non-matched sequences after the search instance(s) performed with BLASTp, DIAMOND or SSEARCH36. |
| hmmscan_[...]_out | HMMER3 output table of the search instance(s) performed against a specific HMMDB. |
| [HMMDB]_non_match.faa | FASTA file with the remaining non-matched sequences after the search instance performed against a specific HMMDB. |
| [...]_omic_table.tbl | Feature table with fields separated by vertical bars. |
| [...]_omic_table.tsv | Feature table with tab-separated fields. |
| arche_report | File which includes the parameters of the run and results. |

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
    
	-v, --verbose	     Use 'yes' to turn on the verbose mode.


# Licence

[GPL v3](https://raw.githubusercontent.com/tseemann/prokka/master/doc/LICENSE.Prokka)

## Author

* Daniel Alonso
* email: gundizalvus16@hotmail.com
