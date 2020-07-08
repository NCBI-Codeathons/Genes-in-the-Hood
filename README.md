# Identify antiphage defense systems in the bacterial pangenome
NCBI Datasets Codeathon Team 3

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
  * for each protein and hmm pair, run `rpstblastn` (B.i) and `tblastn` (B.ii)

* retrieve (protein.faa, genomic.fna, gff) for 80k g-proteobacteria assemblies in Datasets

* extract cas9 (WP) sequences from protein.faa -- single fasta

* extract list of WP accessions

* retrieve TIGRFAM ids associated with WPs, allowing for retrieval of [HMM models](https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR01865.1.HMM)

* hmmer search using cas9 HMM profiles vs genomic nucleotide sequence

