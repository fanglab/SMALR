from setuptools import setup, Extension, find_packages
import sys

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
		print "SMALR requires Python 2.7"
		sys.exit(-1)

globals = {}
execfile("smalr/__init__.py", globals)
__version__ = globals["__version__"]
__author__ = globals["__author__"]
__email__ = globals["__email__"]

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
	'pbcore >= 0.9.4',
	'h5py >= 2.0.1',
	'numpy >= 1.7.1, < 1.14',
	'cython >= 0.17',
	]

setup_requirements = [
	# TODO(jbeaulaurier): put setup requirements (distutils extensions, etc.) here
	"nose",
]

setup(
		name = 'smalr',
		version=__version__,
		author=__author__,
		author_email=__email__,
		description='A pipeline for single molecule-level modification analysis of SMRT reads',
		long_description=readme,
		url='https://github.com/jbeaulaurier/SMALR',
		
		packages=find_packages(include=['smalr']),
		include_package_data=True,
		package_data= {'smalr': ['R/*.r']},
		zip_safe= False,
		keywords='smalr methylation single-molecule bacteria SMRT pacbio',

		entry_points={"console_scripts" : ["smalr = smalr.smalr_multicontig:main"]},
		install_requires=requirements,

		classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Developers',
		'License :: OSI Approved :: BSD License',
		'Natural Language :: English',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
    	],
    	test_suite='nose.collector',
		tests_require=test_requirements,
		setup_requires=setup_requirements,
)
