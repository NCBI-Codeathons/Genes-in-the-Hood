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
(i) create a BLAST db of cas9 proteins
  * using downloaded protein.faa, extract proteins ascribed as "cas9"

(ii) blastp all proteins from "without cas9" assemblies in g-proteobacteria against the cas9 BLAST db &mdash; to find unannotated cas9-like proteins
  * select a representative set of cas9 protein sequences -- to build a query set

(iii) build nucleotide blast database with genomic sequences from all assemblies in g-proteobacteria
  * execute tblastn with the cas9 protein query set

(iv) build protein blast database with a proteins from assemblies in g-proteobacteria
  * execute blastp with the cas9 query set

(v) create a protein profile for known cas9 proteins, use hmmer?
  * domain-based search using rpsblast
    ** look into use of tigrfam as input to makeprofiledb


### (C) Neighborhood analysis
- for genes upstream / downstream of Cas9 (from A ii), are there any patterns?
- Other Cas-related genes?
- For non Cas genes, characterize:  domain structure, GO terms?


## Work Plan 
* retrieve (protein.faa, genomic.fna, gff) for 80k g-proteobacteria assemblies in Datasets
* extract cas9 (WP) sequences from protein.faa -- single fasta
* for B.iii, generate 
