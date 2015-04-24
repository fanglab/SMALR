The SMALR package conducts single-molecule level interrogation of the methylation status of SMRT reads. There are two protocols available for use within the pipeline, SMsn and SMp.

##########################################
SMsn: Single-molecule, single nucleotide analysis
Each motif site on each sequencing molecule is assessed 
for methylation status. This is designed for use with 
short (~250bp) sequencing library preps, where the long 
read lengths of SMRT reads enables multiple passes over 
each motif site. The reliability of the SMsn scores 
increases with more passes (i.e. higher single-molecule 
coverage).
###########################################


###########################################
SMp: Single-molecule, motif-pooled analysis
All motif sites on a sequencing molecule are pooled 
together and the molecule-wide methylation status for 
the given motif is assessed. This is designed for use 
with long (10Kb+) sequencing library preps, where each 
single long subread can span many distinct motif sites. 
The reliability of the SMp scores increases with 
increasing number of distinct motif sites contained in 
the subread.
###########################################

############
Installation:
############
1. Create and activate a python virtual environment with virtualenv
2. Install package dependencies
   - numpy >= 1.9.1
   - h5py >= 2.4.0
   - pbcore >= 0.9.4
   - pysam >= 0.8.1
   - rpy2 >= 2.5.6
   - ???
3. Clone SMALR source code from GitHub
4. Install SMALR inside your virtual environment
   - cd  <SMALR repository>
   - python setup.py install
5. Confirm installation by testing both SMsn and SMp protocols 
   - cd test
   - sh run_test_SMsn.sh
   - sh run_test_SMp.sh

Pipeline output
