SMALR documentation
===================
SMALR: a framework for single-molecule level interrogation of the methylation status of SMRT reads.

The SMALR package conducts single-molecule level interrogation of the methylation status of SMRT reads. There are two protocols available for use within the pipeline, SMsn and SMp.

SMsn: single-molecule, single nucleotide analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Each motif site on each sequencing molecule is assessed for methylation status. This is designed for use with short (~250bp) sequencing library preps, where the long read lengths of SMRT reads enables multiple passes over each motif site. The reliability of the SMsn scores increases with more passes (i.e. higher single-molecule coverage).

SMp: single-molecule, motif-pooled analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All motif sites on a sequencing molecule are pooled together and the molecule-wide methylation status for the given motif is assessed. This is designed for use with long (10Kb+) sequencing library preps, where each single long subread can span many distinct motif sites. The reliability of the SMp scores increases with increasing number of distinct motif sites contained in the subread.


Contribute
^^^^^^^^^^
* Issue tracker: `GitHub <https://github.com/fanglab/smalr/issues>`__
* Source code: `GitHub <https://github.com/fanglab/smalr>`__


Contents
^^^^^^^^
.. toctree::
   :maxdepth: 2

   readme
   installation
   usage
   contributing

Search
^^^^^^

* :ref:`search`