# CRISPRDataCraft_v0.1

  CRISPRDataCraft is a software pipeline designed to enable the identification and statistical analysis of genome editing events. 
  
  CRISPRDataCraft workflow is structured as follows:

    1. Barcode-based Read Splitting: The next-generation sequencing reads are aplit based on barcode sequences provide in ./example/barcode.csv
    2. Paired-end Read Merging: Paired-end reads from next-generation sequencing experiments are merged using the flash2 software (available at https://ccb.jhu.edu/software/FLASH/)
    3. Target Sequence Extraction: The target sequence, flanked by two predefined fragments, is extracted. The corresponding reference sequence and flanking fragments are provided in ./example/gene.fa
    4. Customized Sequence Alignment: The target segments in the sequencing reads are aligned to the reference sequence using an customized algorithm. This alignment algorithm employs a "divide and conquer" strategy: it identifies matching seeds, extends them bidirectionally until mismatches are enountered, and recursively searches for new seeds in the remaining regions, ensuring both efficiency and accuracy in the alignment process.
    5. Mutation Quantification: Mutations, insertions and deletions are quantified to determine whether a read is modified or unmodified by genome editing. This module was developed by modifying CRISPResso2 software (https://github.com/pinellolab/CRISPResso2/tree/master) 
  
  In addition, CRISPRDataCraft provides the following tools, with output of IGDB_split_merge_seed_aligner_local.py as input files:
  
    1. Codon Shift Analysis: Statistical evaluation of codon shifts resulting from genome editing.
       Note: The reference fragment in gene.fa must be a CDS starting at the first position of the protein-encoding triplet codon. 
    2. Desired Mutant Analysis: Statistical assessment of desired mutants generated through genome editing
    
  These functionalities enhance the tool's capability to provide comprehensive insights into genome editing outcomes

# installation
  Compile cython:
  python3 setup_compare.py build_ext --inplace

# Usage
  CRISPRDataCraft:
  
  python3 IGDB_split_merge_seed_aligner_local.py -r1 ./example/1_clean_r1.fq.gz -r2 ./example/1_clean_r2.fq.gz -b ./example/barcode_1.csv -a ./example/gene_1.fa -p 5 -o output_exp1
  
  python3 IGDB_split_merge_seed_aligner_local.py -r1 ./example/2_clean_r1.fq.gz -r2 ./example/2_clean_r2.fq.gz -b ./example/barcode_2.csv -a ./example/gene_2.fa -p 5 -o output_exp2


  Codon Shift Analysis:
  
  python IGDB_analyze_codon_shift.py -pa output_exp1 -o1 ./output_exp1/CSAout1.csv -o2 ./output_exp1/CSAout2.csv  


  Desired Mutant Analysis:
  
  python IGDB_take_desired_mutants.py -aft ./output_exp2/ -dsf ./example/ListbyProduct.txt -o ./output_exp2/Desired_mutants.csv -top 10 -pct 1

