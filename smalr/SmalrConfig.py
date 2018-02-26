import os,sys
import optparse
import logging
import subprocess

class RunnerConfig:
	def __init__( self ):
		self.__parseArgs( )

	def __parseArgs( self ):
		"""Handle command line argument parsing"""

		usage = """%prog [--help] [options] <inputs_file>

		This program will take a native aligned_reads.cmph5 and compare each molecule's
		kinetics profile with a matching WGA aligned_reads.cmp.h5. There are two protocols
		available for use within the pipeline, SMsn and SMp. The single argument for both 
		protocols is an inputs_file containing paths to the relevant files.

		The inputs_file has the following format:
		
		native_cmph5 : /path/to/native/aligned_reads.cmp.h5
		fastq        : /path/to/native/CCS_reads.fastq       (DEPRECATED; specify NONE)
		wga_cmph5    : /path/to/WGA/aligned_reads.cmp.h5     
		ref          : /path/to/sample/reference.fasta

		With the updated support for aligned BAM files (in addition to *.cmp.h5 files), 
		please supply the path to the aligned BAM file using the 'native_cmph5' and 
		'wga_cmph5' lines in the input file.
		

		SMsn: Single-molecule, single nucleotide analysis
		Each motif site on each sequencing molecule is assessed for methylation
		status. This is designed for use with short (~250bp) sequencing library
		preps, where the long read lengths of SMRT reads enables multiple passes
		over each motif site. The reliability of the SMsn scores increases with
		more passes (i.e. higher single-molecule coverage).

		Example usage for SMsn analysis:
		smalr -i --SMsn --motif=GATC --mod_pos=2 --procs=4 -c 5 input.txt


		SMp: Single-molecule, motif-pooled analysis
		All motif sites on a sequencing molecule are pooled together and the
		molecule-wide methylation status for the given motif is assessed. This
		is designed for use with long (10Kb+) sequencing library preps, where each
		single long subread can span many distinct motif sites. The reliability of
		the SMp scores increases with increasing number of distinct motif sites
		contained in the subread.

		Example usage for SMp analysis:
		smalr -i --SMp --motif=GATC --mod_pos=2 --procs=4 -c 5 input.txt
		"""

		parser = optparse.OptionParser( usage=usage, description=__doc__ )

		parser.add_option( "-d", "--debug", action="store_true", help="Increase verbosity of logging" )
		parser.add_option( "-i", "--info", action="store_true", help="Add basic logging" )
		parser.add_option( "--logFile", type="str", help="Write logging to file [log.out]" )
		parser.add_option( "--out", type="str", help="Filename to output SMsn/SMp results [<SMsn/SMp>.out]" )
		parser.add_option( "-c", "--nativeCovThresh", type="int", help="Per mol/strand coverage threshold below which to ignore molecules [10]" )
		parser.add_option( "-m", "--motif", type="str", help="(Required) The sequence motif to be analyzed [None]" )
		parser.add_option( "-s", "--mod_pos", type="int", help="(Required) The modified position in the motif to be analyzed (e.g. for Gm6ATC, mod_pos=2) [None]" )
		parser.add_option( "--wgaCovThresh", type="int", help="Aggregate WGA coverage threshold below which to skip analysis at that position [10]" )
		parser.add_option( "--SMsn", action="store_true", help="Use short-library, single-nucleotide detection protocol. [False]" )
		parser.add_option( "--SMp", action="store_true", help="Use long-library epigenetic phasing protocol (pool IPDs from each subread). [False]" )
		parser.add_option( "--procs", type="int", help="Number of processors to use [4]" )
		parser.add_option( "--align", action="store_true", help="[DEPRECATED; not supported for BAM input] Align native reads to reference to avoid real SNP positions. Only use when expecting sequence heterogeneity in sample (i.e. mtDNA). [False]" )
		parser.add_option( "--upstreamSkip", type="int", help="Number of bases 5' of a CCS-detected, molecule-level SNP to skip in analysis (only when using --align) [10]" )
		parser.add_option( "--downstreamSkip", type="int", help="Number of bases 3' of a CCS-detected, molecule-level SNP to skip in analysis (only when using --align) [10]" )
		parser.add_option( "--minSubreadLen", type="int", help="Minimum length of a subread to analyze [100]" )
		parser.add_option( "--minAcc", type="float", help="Minimum accuracy of a subread to analyze [0.8]" )
		# parser.add_option( "--minMapQV", type="int", help="Minimum mapQV of a subread to analyze [240]" ) #DEPRECATED
		parser.add_option( "--natProcs", type="int", help="Number of processors to use for native molecule analysis [same as procs]" )
		parser.add_option( "--leftAnchor", type="int", help="Number of left bp to exclude around subread-level alignment errors [1]" )
		parser.add_option( "--rightAnchor", type="int", help="Number of right bp to exclude around subread-level alignment errors [1]" )
		parser.add_option( "--write_vars", type="str", help="Write mol-specific variant calls to this filename (only when using --align) [None]" )
		parser.add_option( "--useZMW", action="store_true", help="Index molecules using ZMW/movie ID's instead of molecule IDs (if all alignments all have unique molecule IDs) [False]" )
		parser.set_defaults( logFile="log.out",                \
							 debug=False,                      \
							 info=False,                       \
							 out=".out",                       \
							 motif=None,                       \
							 mod_pos=None,                     \
							 wgaCovThresh=10,                  \
							 align=False,			           \
							 upstreamSkip=10,                  \
							 downstreamSkip=10,                \
							 minSubreadLen=100,                \
							 minAcc=0.8,                       \
							 minMapQV=240,                     \
							 procs=4,                          \
							 natProcs=None,                    \
							 SMsn=False,                       \
							 SMp=False,                        \
							 leftAnchor=1,                     \
							 rightAnchor=1,                    \
							 write_vars=None,                  \
							 useZMW=False)

		self.opts, args = parser.parse_args( )
		if len(args) == 0:
			print usage
			sys.exit()


		if self.opts.align:
			raise Exception("the --align option no longer supported. This has no effect unless your data has significant SNP heterogeneity.")

		# Shift the --mod_pos value into Pythonic 0-based indexing
		self.opts.mod_pos -= 1

		# Hard-coded parameters to avoid effect of adapter sequence on polymerase kinetics
		self.opts.firstBasesToSkip = 15
		self.opts.lastBasesToSkip  = 10

		# These options aren't currently compatible with multicontig support
		self.opts.extractAlignmentOnly = False
		self.opts.nat_aligns_flat      = None
		self.opts.wga_aligns_flat      = None

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

		def check_input_fype(expected_ftype, fname):
			if expected_ftype=="alignments":
				valid_exts  = [".cmp.h5", ".bam"]
				ext_matches = [fname.endswith(x) for x in valid_exts]
				if not any(ext_matches):
					raise Exception("Looking for *.cmp.h5 or *.bam file, found %s" % fname)
				elif fname.endswith(".cmp.h5"):
					self.opts.align_ftype = "cmp"
				elif fname.endswith(".bam"):
					self.opts.align_ftype = "bam"
			
			if expected_ftype=="fasta":
				valid_exts  = ["fa", "fasta", "fsa"]
				ext_matches = [fname.endswith(x) for x in valid_exts]
				if not any(ext_matches):
					raise Exception("Looking for valid fasta file (*.fasta, *.fa, or *.fsa), found %s" % fname)

		contig_dir = os.getcwd()
		orig_dir   = os.path.dirname(self.input_file)
		os.chdir(orig_dir)
		for line in open(self.input_file).xreadlines():
			if line.split(":")[0].strip() == "fastq":
				self.fastq        = os.path.abspath(line.split(":")[1].strip())
			
			elif line.split(":")[0].strip() == "ref":
				self.ref          = os.path.abspath(line.split(":")[1].strip())
				self.opts.ref     = self.ref
				check_input_fype("fasta", self.ref)
			
			elif line.split(":")[0].strip() == "wga_cmph5":
				self.wga_cmph5    = os.path.abspath(line.split(":")[1].strip())
				check_input_fype("alignments", self.wga_cmph5)
			
			elif line.split(":")[0].strip() == "native_cmph5":
				self.native_cmph5 = os.path.abspath(line.split(":")[1].strip())
				check_input_fype("alignments", self.native_cmph5)
			
			else:
				raise Exception("Unexpected field in the input file!\n%s" % line)

		if self.opts.align_ftype=="bam":
			self.opts.useZMW = True

		if not os.path.exists("%s.fai" % self.opts.ref):
			faidx_CMD = "samtools faidx %s" % self.opts.ref
			p         = subprocess.Popen(faidx_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdOutErr = p.communicate()
			sts       = p.returncode
			if sts != 0:
				for entry in stdOutErr:
					print entry
				raise Exception("Failed command: %s" % CMD)

		os.chdir(contig_dir)
