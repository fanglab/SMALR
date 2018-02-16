import os,sys,shutil
import optparse
import logging
from pbcore.io.align.CmpH5IO import *
from pbcore.io import openIndexedAlignmentFile
from itertools import groupby
from smalr import SmalrRunner
import SmalrConfig
from itertools import groupby

class Smalr_multicontig_runner:
	def __init__ ( self ):
		"""

		"""
		self.Config  = SmalrConfig.RunnerConfig( )

	def get_reference_contigs( self, cmph5 ):
		"""
		Pull out the list of contigs in the h5 file.
		"""
		# reader  = CmpH5Reader(cmph5)
		reader  = openIndexedAlignmentFile(cmph5, self.Config.ref)
		contigs = set(map(lambda x: (x[3], x[2]), reader.referenceInfoTable))
		return contigs

	def fasta_iter(self, fasta_name):
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

	def check_contig_names( self, contigs ):
		# Check cmp.h5 contig names against names in the native cmp.h5
		cmp_contig_names = set(map(lambda x: x[0], contigs))
		
		fa_iter          = self.fasta_iter(self.Config.ref)
		#f = open(self.Config.ref+".tmp", "w")
		for name,seq in fa_iter:
			#if name.find("|quiver")>-1:
			#	name = name.split("|")[0]
			#	f.write(">%s\n" % name)
			#	f.write("%s\n" % seq)	
			if name not in cmp_contig_names:
				raise Exception("%s in %s not found in %s!" % (name, self.Config.ref, self.Config.native_cmph5))
			else:
				pass
		#f.close()
		#self.Config.ref = self.Config.ref+".tmp"

	def run( self ):
		nat_contigs  = self.get_reference_contigs(self.Config.native_cmph5)
		self.check_contig_names( nat_contigs )
		
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
			dir_name = "%s_%s" % (nat_contig[0], protocol)
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
