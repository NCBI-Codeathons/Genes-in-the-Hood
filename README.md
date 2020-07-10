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






## The Project

### (A) Retrieve g-proteobacteria assemblies
(i)  retrieve assemblies with GFF, genomic and protein sequences

(ii) find cas9 genes and the gene neighborhood (10 genes upstream/downstream)

(iii) identify/classify genome as "having cas9" vs "without cas9" (based on the symbol "cas9")


### (B) Use cas9 protein/HMM profile sequence to look for similar proteins but without the cas9 symbol
(i) rpstblastn using cas9 HMM profiles (TIGR03031 and TIGR01865)

* build nucleotide blast database with genomic sequences from all assemblies in g-proteobacteria

* extract cas9 protein sequences into a single `cas9.faa`

* retrieve TIGRFAM profiles, build BLAST profile db

* execute rpstblastn with the genomic g-proteobacteria sequences


(ii) tblastn using cas9 proteins

* build a cas9 protein query set &mdash; a subset of `cas9.faa`

* execute tblastn with the `cas9_filtered.faa` protein query set

* create a protein profile for known cas9 proteins, use hmmer?
  * domain-based search using rpsblast
    * look into use of tigrfam as input to makeprofiledb


### (C) Neighborhood analysis

* for genes upstream / downstream of cas9 (from A ii), are there any patterns?

* For non-cas genes, characterize:  domain structure, GO terms?

* cas genes clustered in the neighborhood - Cas1, Cas2, Cas4

* Other cas-related genes?

* neighborhood analysis for test set (tblastn) done - Cas9 is in a cluster with Cas1, Cas2 and Cas4 

* 10 genes on either side vs. 10 kb

## Work Plan 
* experiment: testset-1
  for B, create a small test dataset to test i and ii. small number of assemblies and a selected set of cas9 proteins with associated hmms
  * choose a handful of assemblies -- random selection: GCF_011754595.1 GCF_902706585.1 GCF_011683955.1 GCF_009811535.1 GCF_000948985.1 GCF_002215215.1 GCF_010592905.1
  * build a BLAST genomic nucleotide db -- done
  * extract associated cas9 protein.faa (to `cas9.faa`) -- Ray
  * for each protein accession in `cas9.faa`, download HMM profiles (done)
    * find associated HMM ids, and download the HMM profiles
    * determine NCBI ftp URL for HMM profile
    * download the HMM profiles
  * for each protein and hmm pair, run `rpstblastn` (B.i) and `tblastn` (B.ii) (testing done - for single assembly in B.i; test_cas9.aa input for B.ii)

* retrieve (protein.faa, genomic.fna, gff) for 80k g-proteobacteria assemblies in Datasets

* extract cas9 (WP) sequences from protein.faa -- single fasta 

* extract list of WP accessions

* retrieve the TIGRFAM HMM profiles associated with `cas9` from [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/).
  * [TIGR03031](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR03031.1.HMM) - type II-B CRISPR-associated RNA-guided endonuclease Cas9/Csx12
  * [TIGR01865](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR01865.1.HMM) - type II CRISPR RNA-guided endonuclease Cas9

* hmmer search using cas9 HMM profiles vs genomic nucleotide sequence

    `hmmsearch --tblout orf.10_TIGR0_com_HMMsearch TIGR0_combined.HMM orf.10.test.out`

  * input files
    - `TIGR0_combined.HMM` &mdash; the set of HMM profiles
    - `orf.10.test.out` &mdash; the orf generated from orf-finder
   
  * output file
    - `orf.10_TIGR0_com_HMMsearch` &mdash; the tabulated hmmsearch output
   
It also put the alignment on the screen. We found 14 hits from 10 assemblies, same as that form rpsblast. We found at least one cas9 pseudogene ([ORF1616_NC_007880.1:1261533:1259074](https://www.ncbi.nlm.nih.gov/nuccore/NC_007880.1?report=genbank&from=1258761&to=1263689) annotated as frameshifted cas9 pseudogene, no aa sequence in the downloaded assembly protein.faa).







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

#### Methods evaluated on test set:

* tblastn - protein query vs. translated nucleotide database
* rpstblastn - genomic sequence vs. db of HMM profiles
* ORFfinder/HMM profiles - protein query  
  rpsblast  
   hmmsearch - fastest; used this for cas9 discovery 


#### ORFfinder/hmmsearch

@@@@ ADD COMMANDS here @@@@

* For each assembly, find putative ORFs of at least 300 basepairs...

  `ORFfinder -in genomic.fna -g 11 -s 0 -ml 300 -n t -outfmt 0`

* For each downloaded assembly, we used `ORFfinder` to translate all the ORFs (>=100 AAs) and then used `hmmsearch` to find cas9 homologs.

@@@@ ADD code here@@@@

@@@@ ADD image of html table here@@@@
   







### Neighborhood analysis

- 10 genes on either side of cas9

@@@@ ADD methods here @@@@


- Cas proteins cas1, cas2 and cas4 are clustered with cas9   
   cas1 - endonuclease
   cas2 - spacer acquistion in CRISPR array 
   cas4 - exonuclease; capture of new viral DNA sequences  
   
   
- @@@@ any other genes common to the neighboorhood? @@@@



@@@@ ADD image of table here @@@@


@@@@ ADD image of heat map @@@@

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

- FASTA file of cas9 proteins (WP accessions)
- html table with results of `hmmsearch`
- table of gene neighboorhood
- Heat map of the gene neighboorhood


#### Feedback for the Datasets team

- issues with downloading assemblies
- users are interested in GO (Gene Ontology) annotation










      

