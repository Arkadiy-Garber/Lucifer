# RhoGenie

usage: RhoGenie.py [-h] [-bin_dir BIN_DIR] [-bin_ext BIN_EXT] [-outdir OUTDIR]
                   [-out OUT] [--makeplots MAKEPLOTS] [-hmm_dir HMM_DIR] [--R R]
                   
                   -bin_dir Directory of genomes or assemblies. Must be in FASTA amino acid format.
                   
                   -bin_ext Filename extension for the genomes or assemblies in the bin_dir directory
                   
                   -outdir Name of output directory. This directory will be created.
                   
                   -out Basename for output files (optional)
                   
                   --makeplots (y/n) would you like R plots to be generate automatically? (not recommended for servers where 
                   you do not have super user (sudo) permissions.
                   
                   -hmm_dir Directory of HMMs (i.e. the 'HMMs' directory that comes with this tool)
                   
                   --R Directory of R scripts (i.e. the "rscripts" directory that comes with this tool). Not requried if R       
                   plots are not being generated automatically
                   

Quick note to amateurs: If the program is not in your path, then your console should be in the same directory as the RhoGenie.py executable, and you must run it as "./RhoGenie.py -bin_dir..."
