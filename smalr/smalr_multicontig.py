import os,sys,shutil
import optparse
import logging
from pbcore.io.align.CmpH5IO import *
from itertools import groupby
from smalr import SmalrRunner
import SmalrConfig

def fasta_iter(fasta_name):
	"""
	Given a fasta file, yield tuples of (header, sequence).
	"""
	fh = open(fasta_name)
	# ditch the boolean (x[0]) and just keep the header or sequence since
	# we know they alternate.
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		# drop the ">"
		header = header.next()[1:].strip()
		# join all sequence lines to one.
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq

class Smalr_multicontig_runner:
	def __init__ ( self ):
		"""

		"""
		self.Config  = SmalrConfig.RunnerConfig( )

	def get_reference_contigs( self, cmph5 ):
		"""
		Pull out the list of contigs in the h5 file.
		"""
		reader  = CmpH5Reader(cmph5)
		contigs = set()
		for r in reader:
			contigs.add( (r.referenceInfo[3], r.referenceInfo[2]) )
		return contigs

	def run( self ):
		"""

		"""
		nat_contigs  = self.get_reference_contigs(self.Config.native_cmph5)
		
		print "Preparing to iterate over all contigs in %s" % self.Config.native_cmph5
		for nat_contig in nat_contigs:
			print "	- %s (%s)" % (nat_contig[1], nat_contig[0])

		if self.Config.opts.SMsn:
			protocol = "SMsn"
		elif self.Config.opts.SMp:
			protocol = "SMp"
		else:
			raise Exception("Unknown pipeline protocol! Use --SMsn or --SMp!")

		abs_input_fn = os.path.abspath(self.Config.input_file)
		for i,nat_contig in enumerate(nat_contigs):
			dir_name = "%s_%s" % (nat_contig[1], protocol)
			if os.path.exists(dir_name):
				shutil.rmtree(dir_name)
			os.mkdir(dir_name)
			os.chdir(dir_name)
			runner = SmalrRunner( i, nat_contig, abs_input_fn, self.Config )
			runner.run()
			os.chdir("../")

def main():
	app = Smalr_multicontig_runner()
	sys.exit( app.run() )