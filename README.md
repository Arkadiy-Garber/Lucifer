# Lucifer
Identification of light-sensing and light-producing genes in various 'omics datasets

## Sample command:
    ./Lucifer.py -bin_dir DirectoryWithFASTAs/ -bin_ext faa -outdir OutputDirectoryName/ -out OutputPrefix/ -hmm_dir HMMs/

## Sample command with mapping information provided (via sorted BAM file)
    ./Lucifer.py -bin_dir DirectoryWithFASTAs/ -bin_ext faa -outdir OutputDirectoryName/ -out OutputPrefix/ -hmm_dir HMMs/ -bams MappingInfo.sorted.bam
