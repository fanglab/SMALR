=============
SMALR summary
=============

SMALR: a framework for single-molecule level interrogation of the methylation status of SMRT reads.

The SMALR package conducts single-molecule level interrogation of the methylation status of SMRT reads. There are two protocols available for use within the pipeline, SMsn and SMp.

SMsn: single-molecule, single nucleotide analysis
-------------------------------------------------
Each motif site on each sequencing molecule is assessed for methylation status. This is designed for use with short (~250bp) sequencing library preps, where the long read lengths of SMRT reads enables multiple passes over each motif site. The reliability of the SMsn scores increases with more passes (i.e. higher single-molecule coverage).

SMp: single-molecule, motif-pooled analysis
-------------------------------------------
All motif sites on a sequencing molecule are pooled together and the molecule-wide methylation status for the given motif is assessed. This is designed for use with long (10Kb+) sequencing library preps, where each single long subread can span many distinct motif sites. The reliability of the SMp scores increases with increasing number of distinct motif sites contained in the subread.

Documentation
-------------
For a comprehensive guide on how to install and run mBin, please see the full `documentation <https://smalr.readthedocs.io/en/latest/>`__.

Citation
--------
Beaulaurier J, Zhang XS, Zhu S, Sebra R, Rosenbluh C, Deikus G, Shen N, Munera D, Waldor MK, A Chess, Blaser MJ, Schadt EE, Fang G. `Single molecule-level detection and long read-based phasing of epigenetic variations in bacterial methylomes <http://www.nature.com/articles/ncomms8438>`__. *Nature Communications* **6**, 7438 (2015).