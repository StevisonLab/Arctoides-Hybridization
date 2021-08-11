# dN/dS Analysis

## Extracting coding sequences

The shell script ExtractCodingSequences.sh does most of the heavy lifting here. Details are explained in comments in the script but briefly it requires a VCF file with all samples of interest and a reference FASTA genome. It also requires several third-party programs to run (module loads in shell script) and requires the ReplaceFastaNames R and Python scripts to be present in the same directory. Other than what is in the script, it is worth noting that the VCF file should contain both invariant sites and SNP's, but not indels. We used bcftools view to retain only "ref" and "snp" types (http://samtools.github.io/bcftools/bcftools.html#expressions). The invariant sites and SNP's are both needed to create an accurate variant consensus FASTA sequence. The shell script will mask missing sites but not invariant, so the former should be explicitly given in the VCF file. Indels should be removed because they will throw off alignment given by gffread in the shell script.

## dN/dS

The codeml.sh script can parallelize codeml runs through the gene_sequences folder created from the above processes. This script requires paml to be installed on the computer used.

Details of dNdS_extraction.sh are provided in the script itself. Briefly, it will extract dN/dS values from PAML mlc files (generated from codeml.sh) for each branch (Arc,Sin and Fas) and average their values by topology weight and then average all transcripts for a given gene to generate per-gene estimates of dN/dS. This script requires only commands typical on unix/linux computers (such as awk).