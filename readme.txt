Organohalide-cycling gene database (HaloCycDB)
As the Earth's largest reservoir of dissolved organic carbon, sulfur, and halogen species, the ocean is a hotspot for biogeochemical element cycling, potentially influencing global food webs and climate. Nonetheless, in contrast to the extensively studied organic carbon and sulfur cycles, investigation on the biogeochemical cycling of organohalides in the ocean is still in its infancy. Metagenome sequencing has offered new opportunities to advance our understanding of microbially mediated organohalide cycling. However, there are problems in applying publicly available orthology databases to identify organohalide-cycling genes such as inefficient database searching, inaccurate functional annotation and low coverage of organohalide-cycling genes/pathways. Therefore, it is essential to develop a comprehensive and accurate database for characterizing organohalide-cycling genes. To solve those problems, we constructed a manually curated organohalide-cycling gene database (HaloCycDB) for metagenomic analysis of microbially mediated organohalide-cycling process in the environment. Here, protein sequences from HaloCycDB were recruited from multiple public databases including UniProt, NCBI nr, KEGG, COG, eggNOG, arCOG. HaloCycDB covers 7 dehalogenation processes (oxidative dehalogenase, hydrolytic dehalogenase, reductive dehalogenase, glutathione S-transferase, methyltransferase, dehydrochlorinase, and halohydrin dehalogenase) and 4 halogenation processes (flavin-dependent halogenase, vanadium-dependent halogenase, nonheme iron-dependent halogenase, and S-adenosyl-L-methionine-dependent halogenase), containing sequences of 221 functionally-characterized genes and 187,289 representative homologous genes. HaloCycDB and the associated Python scripts will greatly promote the study of microbially mediated organohalide cycling.


Dependencies Tools
1.Linux system (Python scripts for identifying organohalide-cycling genes have been tested on Ubuntu 20.04)
Python > = 3.7
Biopython >=1.78
2.Database searching tools
usearch: https://www.drive5.com/usearch/download.html
diamond: https://github.com/bbuchfink/diamond/releases
blast: https://ftp.ncbi.nlm.nih.gov/blast/executables/
hmmer: https://github.com/EddyRivasLab/hmmer

User guide
Example
Step 1
Search and identify the organohalide-cycling genes from protein sequences, outputting the fasta file annotated as organohalide-cycling genes along with BLAST--formatted annotation results file.
Command: python haloCycDB_searcher_dev.py -f <search> -i <input protein fasta file> -m <blastp|diamond|usearch> -o <output directory> -d <path to HaloCycDB> -e <evalue> -s <similarity> -p <number of threads>

Detailed explanations: 
-f FUNCTION, --function=FUNCTION: search mode: to search and identify the halogenase/dehalogenase genes from protein sequences.
-i INPUT, --input=INPUT: protein sequences file used for annotation derived from bins or contigs after open reading frame prediction.
-m METHOD, --method=METHOD: method selected for annotation: blastp|diamond|usearch; usearch used as default.
-o OUTPUT, --output=OUTPUT: output directory.
-d DATABASE, --database=DATABASE: path to HaloCycDB.
-e EVALUE, --evalue=EVALUE: evalue set for annotation, default=1e-5.
-s SIMILARITY, --similarity=SIMILARITY: similarity(identity) between query and object sequence, range (0~1), default=0.4.
-p PROCESS, --process=PROCESS: number of processes for annotation, default=4.

Step 2
Map the sequence IDs to organohalide-cycling gene information including host taxa of organohalide-cycling genes, abundance of organohalide-cycling host organisms, halogenation/dehalogenation processes of organohalide-cycling genes, and their associated catalyzing reactions and substrates/products. The Python script halocyc_parse.py is used to generate organohalide-cycling gene information from BLAST-like results against the HaloCycDB database in step 1. 
Note: the abuncance/coverage file calculated using Bowtie2, CheckM, Salmon or other tools and the taxon annotation file for the classification of contigs/bins should be prepared.
Command: python halocyc_parse.py -i <halocyc_genename.txt> -d <path to HaloCycDB>/halocyc.gene -s <bin|contig> -t <taxon annotation file>  -a <abundance file> -o <output file>

Detailed explanations:
-i INPUT, --input=INPUT: input file (halocyc_genename.txt; output file from step1).
-d DATABASE, --database=DATABASE: Database file halocyc.gene.
-s STYLE, --style=STYLE: style, 'bin' or 'contig', consistent with the type of protein sequences file used for annotation in sep1.
-t TAXON, --taxon=TAXON: taxon annotation file (style for contig) or taxon output directory of gtdbtk (style for bin).
-a ABUND, --abund=ABUND: abundance file of bin or contig.
-o OUTPUT, --output=OUTPUT: output file.

Step 3
Merge the annotation result file and sequence file obtained from step 1.
Command: cat halogenase.fa dehalogenase.fa > halocyc.fa; cat halo.search dehalo.search > halocyc.search

Further Strictly filter based on sequence similarity or conserved motifs using Python script motif_search_ident.py, outputting the filtered annotation result file.
Command: python motif_search_ident.py -i <halocyc.fa> -a < halocyc parse file> -b <halocyc.search> -o <output directory> -t <threshold>

Detailed explanations:
-i INPUT, --input=INPUT input fasta file (halocyc.fa).
-a ANOTATE, --annotate=ANNOTATE: halocyc parse file (output file from step2).
-b BRESULT, --bresult=BRESULT: annotation result file (halocyc.search).
-o OUTPUT, --output=OUTPUT output directory, the output file 'halocyc_genename_parse_filter.txt' and 'halocyc_filter.fa' will be obtained in the output directory.
-t THRESHOLD, --threshold=THRESHOLD: threshold of sequence similarity, range (0~1), default=0.4.

