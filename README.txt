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

############
Installation:
############
1. Create and activate a python virtual environment with virtualenv
   - See http://docs.python-guide.org/en/latest/dev/virtualenvs/
2. Install package dependencies inside your virtual environment
   - numpy >= 1.9.1 ('pip install numpy')
   - cython >= 0.17 ('pip install cython')
   - h5py >= 2.4.0 ('pip install h5py')
   - pbcore >= 0.9.4 (https://github.com/PacificBiosciences/pbcore)
     * Requires Python >= 2.7
3. Clone SMALR source code from GitHub to a local SMALR repository
   - https://github.com/jbeaulaurier/SMALR
4. Install SMALR inside your virtual environment
   - cd  <SMALR repository>
   - python setup.py install
5. Confirm installation by testing both SMsn and SMp protocols 
   - cd test
   - sh run_test_SMsn.sh
     * Should generate the folder ref000001_SMsn, containing pipeline output
   - sh run_test_SMp.sh
     * Should generate the folders ref000001_SMp through ref000004_SMp,
       containing pipeline output

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
only one contig, the results will be placed in ref000001. These results include a 
log detailing the analysis of that contig, the motif positions in that contig (forward 
and reverse strand), a fasta file of that contig, and a results file (SMsn.out or 
SMp.out). This results file contains the following informtation:

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
