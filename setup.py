from setuptools import setup, Extension, find_packages
import sys

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
		print "SMALR requires Python 2.7"
		sys.exit(-1)

globals = {}
execfile("smalr/__init__.py", globals)
__VERSION__ = globals["__VERSION__"]

setup(
		name = 'smalr',
		version=__VERSION__,
		author='John Beaulaurier',
		author_email='john.beaulaurier@mssm.edu',
		description='A pipeline for single molecule-level modification analysis of SMRT reads',
		url='https://github.com/jbeaulaurier/SMALR',
		
		packages= ['smalr'],
		package_dir= {'smalr' : 'smalr'},
		package_data= {'smalr': ['R/*.r']},
		zip_safe= False,
		scripts= ['bin/smalr',
			  'bin/call_smalr.py'],
		install_requires=['pbcore >= 0.9.4',
				  'h5py >= 2.0.1',
				  'numpy >= 1.6.0',
				  'cython >= 0.17'])
