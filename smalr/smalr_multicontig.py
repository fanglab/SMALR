import os,sys,shutil
import optparse
import logging
from pbcore.io.align.CmpH5IO import *
from itertools import groupby
from smalr import SmalrRunner


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
		self.input_file = sys.argv[-1]
		for line in open(self.input_file).xreadlines():
			if line.split(":")[0].strip() == "native_cmph5":
				self.native_cmph5 = line.split(":")[1].strip()

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
		nat_contigs  = self.get_reference_contigs(self.native_cmph5)
		
		abs_input_fn = os.path.abspath(self.input_file)
		for i,nat_contig in enumerate(nat_contigs):
			if os.path.exists(nat_contig[1]):
				shutil.rmtree(nat_contig[1])
			os.mkdir(nat_contig[1])
			os.chdir(nat_contig[1])
			runner = SmalrRunner( i, nat_contig, abs_input_fn )
			runner.run()
			os.chdir("../")

# def main():
if __name__ == "__main__":
	app = Smalr_multicontig_runner()
	sys.exit( app.run() )