# Antiphage CRISPR-Cas defense system in gamma-protobacterial genomes: Identify Cas9 proteins and their gene neighboorhood 
# Background
Bacteriophages (phages) are viruses that infect bacteria. There are multiple antiphage systems employed by bacteria to protect themselves from phages. One of these is the CRISPR-Cas system, which is a type of acquired immunity. The two components of this system are:
- CRISPR (**C**lustered **R**egularly **I**nterspaced **S**hort **P**alindromic **R**epeats)

   This is a specialized region of the bacterial genome that contains repeats and spacers. The spacers consist of virus sequences from previous phage infections that are used to recognize future attacks 
- Cas (**C**RISPR **ass**ociated) proteins  

   Cas proteins use CRISPR sequences as guides to recognize and cleave specific strands of virus DNA. Cas9 is an endonuclease (cleaves internal bonds in DNA) that causes site-directed double-stranded breaks.    

The CRISPR-Cas9 system also has applications in basic research, and in biotechnology, including genome editing for treatment of genetic disorders.  

# Hypothetical user story 
A research group would like to identify Cas9 and Cas9-like proteins in gammaproteobacterial genomes, and download genes in the neighborhood of Cas9. They start their search with the gene symbol “Cas9” with the goal of finding and downloading
- protein FASTA files for annotated Cas9 proteins (e.g., WP_ proteins in RefSeq) in gammaproteobacteria (taxonomic rank = Class; txid 1236; NCBI BLAST name: g-proteobacteria)
- protein FASTA files for unannotated Cas9-like proteins (similar proteins without symbols or informative names)
- protein FASTA for genes in the neighborhood of Cas9 (10kb on each side)
- GO (Gene Ontology) terms associated with neighboring genes  

How many of the neighborhoods are unique? Can they be classified and made non-redundant?


## The Team
- Vinita Joardar
- Peter Meric
- Jianli Dai
- Greg Schuler
- Wayne Matten
- Ray Anderson


## The Project


### (A) Retrieve g-proteobacteria assemblies
(i)  retrieve assemblies with GFF, genomic and protein sequences

(ii) find cas9 genes and the gene neighborhood (10kB upstream/downstream)

(iii) identify/classify genome as "having cas9" vs "without cas9" (based on the symbol "cas9")


### (B) Use cas9 protein sequence to look for similar proteins but without the cas9 symbol
(i) tblastn using cas9 HMM profiles

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
- for genes upstream / downstream of Cas9 (from A ii), are there any patterns?
- Other Cas-related genes?
- For non Cas genes, characterize:  domain structure, GO terms?


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

* retrieve TIGRFAM ids associated with WPs, allowing for retrieval of [HMM models](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR01865.1.HMM) (done)

* hmmer search using cas9 HMM profiles vs genomic nucleotide sequence

* neighborhood analysis for test set (tblastn) done - Cas9 is in a cluster with Cas1, Cas2 and Cas4 


## Results
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

### Identifying unannotated Cas9 

Methods evaluated on test set:

* rpstblastn

* tblastn

* ORFfinder/rpsblast

* ORFfinder/hmm-search - fastest 

### Neighborhood analysis

* Cas proteins present - Cas1, Cas2, Cas4

* 10 genes on either side vs. 10 kb
