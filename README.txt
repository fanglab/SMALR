The SMALR package conducts single-molecule level interrogation of the 
methylation status of SMRT reads. There are two protocols available 
for use within the pipeline, SMsn and SMp.

SMsn: Single-molecule, single nucleotide analysis
Each motif site on each sequencing molecule is assessed for methylation 
status. This is designed for use with short (~250bp) sequencing library 
preps, where the long read lengths of SMRT reads enables multiple passes 
over each motif site. The reliability of the SMsn scores increases with 
more passes (i.e. higher single-molecule coverage).

SMp: Single-molecule, motif-pooled analysis
All motif sites on a sequencing molecule are pooled together and the 
molecule-wide methylation status for the given motif is assessed. This 
is designed for use with long (10Kb+) sequencing library preps, where each 
single long subread can span many distinct motif sites. The reliability of 
the SMp scores increases with increasing number of distinct motif sites 
contained in the subread.

###########
Dependencies
###########
- python >= 2.7
- bwa >= 0.7
- samtools >= 0.1.18
- R >= 3.0

############
Installation:
############
1. Create and activate a python virtual environment using virtualenv
   - Detailed instructions at http://docs.python-guide.org/en/latest/dev/virtualenvs/
   - Two routes
      * If pip is installed:
         1) pip install virtualenv
         2) virtualenv smalr_venv
         3) . smalr_venv/bin/activate
      * If pip is NOT installed:
         1) curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-13.0.1.tar.gz
         2) tar xvfz virtualenv-13.0.1.tar.gz 
         3) cd virtualenv-13.0.1
         4) python virtualenv.py smalr_venv
         5) . smalr_venv/bin/activate
   - IMPORTANT: This virtual environment must be activated to install & run SMALR!
2. Clone SMALR source code from GitHub to a local SMALR repository
   - git clone https://github.com/fanglab/SMALR.git
2. Install SMALR and required packages inside your virtual environment
   - cd SMALR
   - ./install.sh
5. Confirm successful installation by testing both SMsn and SMp protocols 
   - cd test
   - ./run_test_SMsn.sh
     * Should generate the contig-specific folder J99_SMsn, containing pipeline output.
   - ./run_test_SMp.sh
     * Should generate the contig-specific folders scf7180000000008|quiver_SMp,
       scf7180000000009|quiver_SMp, scf7180000000010|quiver_SMp, and scf7180000000011|quiver_SMp, 
       each containing pipeline output.

############
Pipeline input
############
Both the SMsn and SMp protocols require an input_files.txt argument:
   - Specifies relevant file paths in the following format:
      native_cmph5: <path to native cmp.h5>
      fastq:        <path to native CCS fastq file> (optional, can specify NONE)
      wga_cmph5:    <path to WGA cmp.h5> (optional, can specify native cmp.h5 if not 
                    available)
      ref:          <path to reference that matches that used in the cmp.h5 files>
   - IMPORTANT NOTES:
      * For good results, cmp.h5 files generated from short library sequencing should 
        be paired with the SMsn protocol.
      * Conversely, cmp.h5 files containing long library-sequenced reads are best when
        using the SMp protocol.
      * A CCS fastq file is only REQUIRED when the --align option is specified. This 
        calls an alignment step that is used to mask out sequencing errors on each 
        molecule. --align should never be used with the SMp protocol and long libraries, 
        as CCS sequencing only works with short libraries.

############
Pipeline output
############
One output directory will be created for each contig in the reference. If there is
only one contig, the results will be placed in the folder named for that contig. These 
results include a log detailing the analysis of that contig, the motif positions in 
that contig (forward and reverse strand), a fasta file of that contig, and a results 
file (SMsn.out or SMp.out). This results file contains the following informtation:

Column  Meaning
1       Contig strand
2       Contig motif position (for SMp, pooled motif sites are summarized by smallest 
        site position)
3       SMsn or SMp score (native score - WGA score)
4       Molecule ID
5       Native score (mean of subread-normalized ln(IPD) values; site- and molecule-
        specific)
6       WGA score (mean of subread-normalized ln(IPD) values; site-specific accross 
        all WGA molecules)
7       Number of data points used to get the molecule-level native score
8       Number of data points used to get the aggregate WGA score
9       Mean length of subreads from the native molecule
