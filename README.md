# Antiphage CRISPR-Cas defense system in gamma-protobacterial genomes: Identify Cas9 proteins and their gene neighborhood 

Vinita Joardar, Peter Meric, Ray Anderson, Jianli Dai, Wayne Matten, Greg Schuler


# Introduction
## Background
Bacteriophages (phages) are viruses that infect bacteria. There are multiple antiphage systems employed by bacteria to protect themselves from phages. One of these is the CRISPR-Cas system, which is a type of acquired immunity. The two components of this system are:
- CRISPR (**C**lustered **R**egularly **I**nterspaced **S**hort **P**alindromic **R**epeats)

   This is a specialized region of the bacterial genome that contains repeats and spacers. The spacers consist of virus sequences from previous phage infections that are used to recognize future attacks 
- Cas (**C**RISPR-**as**sociated) proteins  

Cas proteins use CRISPR sequences as guides to recognize and cleave specific strands of virus DNA. Cas9 is an endonuclease (cleaves internal bonds in DNA) that causes site-directed double-stranded breaks.    

The CRISPR-Cas9 system also has applications in basic research, and in biotechnology, including genome editing for treatment of genetic disorders. 
   
## Taxonomy 

- Taxonomy rank: Class

- Gammaproteobacteria - diverse 

    - human pathogens – E. coli, Salmonella

    - plant pathogens – Pseudomanas, Xanthomonas

    - photoautotrophs – purple sulfur bacteria

    - methane oxidizers – Methylococcus



 

## Hypothetical user story 
A research group would like to identify **Cas9** and Cas9-like proteins in gammaproteobacterial genomes, and download genes in the neighborhood of Cas9. They start their search with the gene symbol “cas9” with the goal of finding and downloading
- protein FASTA files for annotated cas9 genes (e.g., WP_ proteins in RefSeq) in gammaproteobacteria (taxonomic rank = Class; txid 1236; NCBI BLAST name: g-proteobacteria)

- protein FASTA files for unannotated cas9-like proteins (similar proteins without symbols or informative names)

- protein FASTA for genes in the neighborhood of cas9 (10 genes on each side)

- GO (Gene Ontology) terms associated with neighboring genes  

How many of the neighborhoods are unique? Can they be classified and made non-redundant?






## Project Overview

### (A) Assembly download and classification
(i)  retrieve assemblies with GFF, genomic and protein sequences

(ii) find cas9 genes and the gene neighborhood (10 genes upstream/downstream)

(iii) identify/classify genome as "having cas9" vs "without cas9" (based on the symbol "cas9")


### (B) Identify unannotated Cas9 proteins
(i) test a few methods

(ii) settled on running ORFfinder on the genomic.fna files, then passing that to hmmsearch


### (C) Neighborhood analysis

* for genes upstream / downstream of cas9 (from A ii), are there any patterns?

* For non-cas genes, characterize:  domain structure, GO terms?

* cas genes clustered in the neighborhood - Cas1, Cas2, Cas4

* Other cas-related genes?

* neighborhood analysis for test set (tblastn) done - Cas9 is in a cluster with Cas1, Cas2 and Cas4 

* 10 genes on either side vs. 10 kb


## Methods and Results
### Assembly download and classification


The NCBI Datasets command-line tool was used to download all available reference assemblies within the clade for gamaproteobacteria (tax ID 1236).  In total this amounted to 79,802 datasets.                  

An initial classification of these datasets involved scanning the genomic.gff file for the presence or absence of "cas9", annotated either as a gene or a pseudogene.                                           

| category            | assemblies |
|---------------------|---------|
| has cas9 gene       |     418 |
| has cas9 pseudogene |     251 |
| no cas9 annotated   |  79,057 |
| missing genomic.gff |      76 |

For the 76 cases in which the downloaded package did not contain a genomic.gff file, it was noted that the file was present on the FTP site, so this is an apparent bug for NCBI Datasets.                      

In addition, we found two cases which could not be downloaded as fully hydrated packages (after repeated attempts), but could be downloaded in dehydrated form and successfully rehydrated.   

### Identifying unannotated Cas9 proteins

#### Methods evaluated on a test set

* `tblastn` &mdash; protein query vs. translated nucleotide database
* `rpsblast` &mdash; protein sequence vs database of HMM profiles
* `rpstblastn` &mdash; genomic sequence vs. database of HMM profiles
* `ORFfinder` &mdash; to generate putative protein sequences
* `hmmsearch` &mdash; search HMM profile(s) against a sequence database

`hmmsearch` plus `ORFfinder` was deemed the fastest and with sufficient sensitivity. Therefore it was chosen for this `cas9` discovery work.


#### ORFfinder/hmmsearch


* For each assembly, find putative ORFs of at least 300 basepairs

  `ORFfinder -in genomic.fna -g 11 -s 0 -ml 300 -n t -outfmt 0`
  
  where:

   `-g 11` &mdash; The Bacterial, Archaeal and Plant Plastid Code
   
   `-s 0`  &mdash; use ATG only for the ORF start codon
   
   `-n t`  &mdash; Ignore nested ORFs (completely placed within another)

* retrieve the TIGRFAM HMM profiles associated with `cas9` from [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/hmm/current/).
  * [TIGR03031](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR03031.1.HMM) &mdash; type II-B CRISPR-associated RNA-guided endonuclease Cas9/Csx12
  * [TIGR01865](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR01865.1.HMM) &mdash; type II CRISPR RNA-guided endonuclease Cas9

* For each downloaded assembly, we used `ORFfinder` to translate all the ORFs (>=100 AAs) and then used `hmmsearch` to find cas9 homologs.

* hmmer search using cas9 HMM profiles vs genomic nucleotide sequence

    `hmmsearch --tblout orf.10_TIGR0_com_HMMsearch TIGR0_combined.HMM orf.10.test.out`

  * input
    - `TIGR0_combined.HMM` &mdash; the set of HMM profiles
    - protein FASTA of putative ORFs, generated by `ORFfinder`
   
[See the results table](https://htmlpreview.github.io/?https://github.com/NCBI-Codeathons/Genes-in-the-Hood/blob/main/src/html/hmmsearch.results.html)




### Neighborhood analysis

Retrieve the neighborhood of 10 genes/proteins on either side of `cas9`

- Cas proteins `cas1`, `cas2` and `cas4` are clustered with `cas9`
   * `cas1` &mdash; endonuclease
   * `cas2` &mdash; spacer acquistion in CRISPR array 
   * `cas4` &mdash; exonuclease; capture of new viral DNA sequences  


#### Heat maps depicting relative frequency of genes/proteins

[Global gene frequency](https://github.com/NCBI-Codeathons/Genes-in-the-Hood/blob/main/cas9-global-freqs.md)

[Column-based gene frequency](https://github.com/NCBI-Codeathons/Genes-in-the-Hood/blob/main/cas9-column-freqs.md)


#### Top 10 occurring genes/proteins

- cas9
- WP_003016533.1
- WP_003016538.1
- WP_003016535.1
- WP_003016544.1
- WP_003019268.1
- WP_003016518.1
- WP_003016545.1
- WP_003016540.1
- WP_003016537.1


### Conclusions

#### Output for the user

- html table with results of `hmmsearch`
- heat map of the gene neighborhood
- table of gene neighborhood
- _pending_: FASTA file of `cas9` proteins (WP accessions)


#### Feedback for the Datasets team

- issues with downloading assemblies
- users are interested in GO (Gene Ontology) annotation










      

