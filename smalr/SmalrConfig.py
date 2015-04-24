import os,sys
import optparse
import logging

class RunnerConfig:
	def __init__( self ):
		self.__parseArgs( )

	def __parseArgs( self ):
		"""Handle command line argument parsing"""

		usage = """%prog [--help] [options] <inputs_file>

		This program will take a native aligned_reads.cmph5 and compare each molecule's
		kinetics profile with the position/strand-matched control kinetics profile.
		The control kinetics profile is derived from either a sample-matched WGA 
		aligned_reads.cmp.h5 or a context-dependant in silico kinetics profile based
		on a kinetics model that was created using training training data.

		The inputs_file contains full paths to the necessary files in the following format:
		
		fastq        : /path/to/native/CCS_reads.fastq       (optional)
		native_cmph5 : /path/to/native/aligned_reads.cmp.h5
		wga_cmph5    : /path/to/WGA/aligned_reads.cmp.h5     
		ref          : /path/to/sample/reference.fasta
		"""

		parser = optparse.OptionParser( usage=usage, description=__doc__ )

		parser.add_option( "-d", "--debug", action="store_true", help="Increase verbosity of logging" )
		parser.add_option( "-i", "--info", action="store_true", help="Add basic logging" )
		parser.add_option( "-p", "--profile", action="store_true", help="Profile code" )
		parser.add_option( "--logFile", type="str", help="Write logging to file [log.out]" )
		parser.add_option( "--out", type="str", help="Filename to output SMsn/SMp results [<SMsn/SMp>.out]" )
		parser.add_option( "-c", "--nativeCovThresh", type="int", help="Per mol/strand coverage threshold below which to ignore molecules [10]" )
		parser.add_option( "-m", "--motif", type="str", help="(Required) The sequence motif to be analyzed [None]" )
		parser.add_option( "-s", "--mod_pos", type="int", help="(Required) The modified position (0-based) in the motif to be analyzed (e.g. for Gm6ATC, mod_pos=1) [None]" )
		parser.add_option( "--wgaCovThresh", type="int", help="Aggregate WGA coverage threshold below which to skip analysis at that position [10]" )
		parser.add_option( "--align", action="store_true", help="Align native reads to reference to avoid real SNPs positions [False]" )
		parser.add_option( "--firstBasesToSkip", type="int", help="Number of beginning bases in the subread to skip [15]" )
		parser.add_option( "--lastBasesToSkip", type="int", help="Number of final bases in the subread to skip [10]" )
		parser.add_option( "--noCompare", action="store_true", help="Only load and characterize native reads [False]]" )
		parser.add_option( "--nat_aligns_flat", type="str", help="Skip the native cmp.h5 extraction step and use this pre-existing flat file containing alignments instead [None]" )
		parser.add_option( "--wga_aligns_flat", type="str", help="Skip the WGA cmp.h5 extraction step and use this pre-existing flat file containing alignments instead [None]" )
		parser.add_option( "--extractAlignmentOnly", action="store_true", help="Stop the pipeline after extracting the alignments data from the native and WGA cmp.h5 files. [False]" )
		parser.add_option( "--upstreamSkip", type="int", help="Number of bases 5' of a detected SNP/error to skip in analysis [10]" )
		parser.add_option( "--downstreamSkip", type="int", help="Number of bases 3' of a detected SNP/error to skip in analysis [10]" )
		parser.add_option( "--minSubreadLen", type="int", help=" [100]" )
		parser.add_option( "--minAcc", type="float", help=" [0.8]" )
		parser.add_option( "--minMapQV", type="int", help=" [240]" )
		parser.add_option( "--procs", type="int", help="Number of subprocesses to spawn [8]" )
		parser.add_option( "--natProcs", type="int", help="Number of subprocesses to spawn for native molecule analysis [procs]" )
		parser.add_option( "--chunkSize", type="int", help="Number of molecules to process at a time [500]" )
		parser.add_option( "--SMsn", action="store_true", help="For short-library single-nucleotide detection. [False]" )
		parser.add_option( "--SMp", action="store_true", help="For long-library epigenetic phasing (pool IPDs from each subread). [False]" )
		parser.add_option( "--nat_lib", type="str", help="Either 'short' or 'long' [short]" )
		parser.add_option( "--wga_lib", type="str", help="Either 'short' or 'long' [short]" )
		parser.add_option( "--leftAnchor", type="int", help="Left buffer to exclude around mismatches and indels [1]" )
		parser.add_option( "--rightAnchor", type="int", help="Right buffer to exclude around mismatches and indels [1]" )
		parser.add_option( "--write_vars", type="str", help="(Must also use --align) Write mol-specific variant calls to this file. [None]" )
		parser.set_defaults( logFile="log.out",                \
							 debug=False,                      \
							 info=False,                       \
							 profile=False,                    \
							 out=".out",                       \
							 motif=None,                       \
							 mod_pos=None,                     \
							 wgaCovThresh=10,                  \
							 align=False,			           \
							 firstBasesToSkip=15,			   \
							 lastBasesToSkip=10,		       \
							 noCompare=False,                  \
							 nat_aligns_flat=None,             \
							 wga_aligns_flat=None,             \
							 extractAlignmentOnly=False,       \
							 upstreamSkip=10,                  \
							 downstreamSkip=10,                \
							 minSubreadLen=100,                \
							 minAcc=0.8,                       \
							 minMapQV=240,                     \
							 procs=8,                          \
							 natProcs=None,                    \
							 SMsn=False,                       \
							 SMp=False,                        \
							 wga_lib="short",                  \
							 nat_lib="short",                  \
							 leftAnchor=1,                     \
							 rightAnchor=1,                    \
							 write_vars=None,                  \
							 chunkSize=500)

		self.opts, args = parser.parse_args( )

		if (not self.opts.SMsn and not self.opts.SMp) or (self.opts.SMsn and self.opts.SMp):
			raise Exception("Specify SMALR mode using EITHER --SMsn or --SMp (but not both!)")

		if self.opts.SMsn:
			self.opts.out = "SMsn"+self.opts.out
		elif self.opts.SMp:
			self.opts.out = "SMp"+self.opts.out

		if self.opts.motif==None or self.opts.mod_pos==None:
			raise Exception("Please specify both the motif (--motif=<str>) and modified position in the motif (--mod_pos=<int>) to be analyzed!")

		if self.opts.write_vars!=None and not self.opts.align:
			print "WARNING: --align not specified... ignoring --write_vars=%s" % self.opts.write_vars
			self.opts.write_vars = None

		if self.opts.natProcs == None:
			self.opts.natProcs = self.opts.procs

		if len(args) != 1:
			parser.error( "Expected 1 argument." )
		else:
			self.input_file = os.path.abspath(args[0])

		if self.opts.wga_aligns_flat != None and not os.path.exists(self.opts.wga_aligns_flat):
			raise Exception("%s file not found!" % self.opts.wga_aligns_flat)
		if self.opts.nat_aligns_flat != None and not os.path.exists(self.opts.nat_aligns_flat):
			raise Exception("%s file not found!" % self.opts.nat_aligns_flat)

		contig_dir = os.getcwd()
		orig_dir   = os.path.dirname(self.input_file)
		os.chdir(orig_dir)
		for line in open(self.input_file).xreadlines():
			if line.split(":")[0].strip() == "fastq":
				self.fastq        = os.path.abspath(line.split(":")[1].strip())
			elif line.split(":")[0].strip() == "ref":
				self.ref          = os.path.abspath(line.split(":")[1].strip())
			elif line.split(":")[0].strip() == "wga_cmph5":
				self.wga_cmph5    = os.path.abspath(line.split(":")[1].strip())
			elif line.split(":")[0].strip() == "native_cmph5":
				self.native_cmph5 = os.path.abspath(line.split(":")[1].strip())
			else:
				raise Exception("Unexpected field in the input file!\n%s" % line)
		os.chdir(contig_dir)
