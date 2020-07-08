# Identify antiphage defense systems in the bacterial pangenome
NCBI Datasets Codeathon Team 3

## The Team
- Team Member 1
- Team Member 2
- Vinita Joardar
- Peter Meric
- Jianli Dai
- Greg Schuler


## The Project


### (A) Retrieve g-proteobacteria assemblies
(i)  retrieve assemblies with GFF, genomic and protein sequences

(ii) find cas9 genes and the gene neighborhood (10kB upstream/downstream)

(iii) identify/classify genome as "having cas9" vs "without cas9" (based on the symbol "cas9")


### (B) Use cas9 protein sequence to look for similar proteins but without the cas9 symbol
* build nucleotide blast database with genomic sequences from all assemblies in g-proteobacteria

* extract cas9 protein sequences into a single `cas9.faa`

* build a cas9 protein query set &mdash; a subset of `cas9.faa`
  * use tigrfam 
https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR01865.1.HMM

* execute tblastn with the `cas9_filtered.faa` protein query set

* create a protein profile for known cas9 proteins, use hmmer?
  * domain-based search using rpsblast
    * look into use of tigrfam as input to makeprofiledb


### (C) Neightbohood analysis
- for genes upstream / downstream of Cas9 (from A ii), are there any patterns?
- Other Cas-related genes?
- For non Cas genes, characterize:  domain structure, GO terms?


## Work Plan 
* retrieve (protein.faa, genomic.fna, gff) for 80k g-proteobacteria assemblies in Datasets
* extract cas9 (WP) sequences from protein.faa -- single fasta
* extract list of WP accessions
* retrieve TIGRFAM ids associated with WPs, allowing for retrieval of [HMM models](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR01865.1.HMM)
