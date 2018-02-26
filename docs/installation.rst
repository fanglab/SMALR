.. highlight:: shell

============
Installation
============

Fundamental dependencies
------------------------
* python 2.7
* samtools >= 0.1.18
* R >= 3.0

Setting up virtualenv
---------------------
Create and activate a python virtual environment using ``virtualenv``. Detailed instructions `here <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`__.

If pip is installed, ``virtualenv`` can be installed using ``pip``:

.. code-block:: console
	
	$ pip install virtualenv
	$ virtualenv smalr_venv
	$ . smalr_venv/bin/activate

If ``pip`` is NOT installed:
            
.. code-block:: console
	
	$ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-13.0.1.tar.gz
	$ tar xvfz virtualenv-13.0.1.tar.gz 
	$ cd virtualenv-13.0.1
	$ python virtualenv.py smalr_venv
	$ . smalr_venv/bin/activate

**NOTE:** This virtual environment must be activated to install & run SMALR

Installing mBin
---------------

With the virtual environment activated, clone SMALR source code from GitHub to a local SMALR repository:

.. code-block:: console
	
	$ git clone https://github.com/fanglab/SMALR.git

Install SMALR and required packages inside your virtual environment:

.. code-block:: console

	$ cd SMALR
	$ python setup.py install

Confirm successful installation by testing both SMsn and SMp protocols:
   
.. code-block:: console

	$ cd test
	$ ./run_test_SMsn.sh
    
This should generate the contig-specific folder ``J99_SMsn``, containing pipeline output.

.. code-block:: console
	
	$ ./run_test_SMp.sh

This should generate the contig-specific folders ``scf7180000000008|quiver_SMp``, ``scf7180000000009|quiver_SMp``, ``scf7180000000010|quiver_SMp``, and ``scf7180000000011|quiver_SMp``, each containing pipeline output.