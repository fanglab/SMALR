import os,sys,time
import parse_mol_aligns
import subprocess
import optparse
import logging
import numpy as np
from itertools import groupby
from collections import defaultdict, Counter
import math
import pickle
from pbcore.io.align.CmpH5IO import *
import multiprocessing
import glob
import random
import h5py
import rpy2
import copy
import operator
import matplotlib.pyplot as plt

COMP = {"A":"T", "T":"A", "C":"G", "G":"C"}

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

def check_file( fn ):
	if not os.path.exists(fn) or os.path.getsize(fn)==0:
		raise Exception( "%s does not exists or is empty!" % fn )

def run_command( CMD ):
	p         = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdOutErr = p.communicate()
	sts       = p.returncode
	if sts != 0:
		for entry in stdOutErr:
			print entry
		raise Exception("Failed command: %s" % CMD)

def cat_list_of_files( in_fns, out_fn, header=None ):
	"""
	Given a list of filenames, cat them together into a single file. Then cleans up pre-catted
	single files.
	"""
	if len(in_fns)==0:
		raise Exception("There are no test.out chunked files to cat!")
	
	if header != None:
		f = open(out_fn,"w")
		f.write(header)
		f.close()
		cat_CMD  = "cat %s  >> %s" % (" ".join(in_fns), out_fn)
	else:
		cat_CMD  = "cat %s  > %s" % (" ".join(in_fns), out_fn)
	p         = subprocess.Popen(cat_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdOutErr = p.communicate()
	sts       = p.returncode
	if sts != 0:
		raise Exception("Failed cat command: %s" % cat_CMD)
	for fn in in_fns:
		try:
			os.remove(fn)
		except OSError:
			pass
	return out_fn

class Consumer(multiprocessing.Process):
	def __init__(self, task_queue, result_queue):
		multiprocessing.Process.__init__(self)
		self.task_queue   = task_queue
		self.result_queue = result_queue

	def run(self):
		proc_name = self.name
		while True:
			next_task = self.task_queue.get()
			if next_task is None:
				# Poison pill means shutdown
				logging.info("%s: Exiting" % proc_name)
				self.task_queue.task_done()
				break
			logging.info("%s: Starting" % proc_name)
			answer = next_task()
			self.task_queue.task_done()
			self.result_queue.put(answer)
		return

class Smalr_runner:
	def __init__ ( self ):
		"""
		Parse the options and arguments, then instantiate the logger. 
		"""
		self.__parseArgs( )
		self.__initLog( )

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
		parser.add_option( "--contig_name", type="str", help="Name of contig in reference to look at [ref000001]" )
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
							 contig_name="ref000001",          \
							 SMsn=False,                       \
							 SMp=False,                        \
							 wga_lib="short",                  \
							 nat_lib="short",                  \
							 leftAnchor=1,                     \
							 rightAnchor=1,                    \
							 write_vars=False,                 \
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

		if self.opts.wga_aligns_flat != None and not os.path.exists(self.opts.wga_aligns_flat):
			raise Exception("%s file not found!" % self.opts.wga_aligns_flat)
		if self.opts.nat_aligns_flat != None and not os.path.exists(self.opts.nat_aligns_flat):
			raise Exception("%s file not found!" % self.opts.nat_aligns_flat)

		input_file = args[0]
		for line in open(input_file).xreadlines():
			if line.split(":")[0].strip() == "fastq":
				self.fastq        = line.split(":")[1].strip()
			elif line.split(":")[0].strip() == "ref":
				self.ref          = line.split(":")[1].strip()
			elif line.split(":")[0].strip() == "wga_cmph5":
				self.wga_cmph5    = line.split(":")[1].strip()
			elif line.split(":")[0].strip() == "native_cmph5":
				self.native_cmph5 = line.split(":")[1].strip()
			else:
				raise Exception("Unexpected field in the input file!\n%s" % line)

		for i,entry in enumerate(fasta_iter(self.ref)):
			self.ref_size = len(entry[1])
			self.ref_seq  = entry[1]
			# FIXME: need to support multiple contigs!
			break
			if i>0:
				raise Exception("Multi-contig references not supported!")

	def __initLog( self ):
		"""Sets up logging based on command line arguments. Allows for three levels of logging:
		logging.error( ): always emitted
		logging.info( ) : emitted with --info or --debug
		logging.debug( ): only with --debug"""

		if os.path.exists(self.opts.logFile):
			os.remove(self.opts.logFile)

		logLevel = logging.DEBUG if self.opts.debug \
					else logging.INFO if self.opts.info \
					else logging.ERROR

		self.logger = logging.getLogger("")
		self.logger.setLevel(logLevel)
		
		# create file handler which logs even debug messages
		fh = logging.FileHandler(self.opts.logFile)
		fh.setLevel(logLevel)
		
		# create console handler with a higher log level
		ch = logging.StreamHandler()
		ch.setLevel(logLevel)
		
		# create formatter and add it to the handlers
		logFormat = "%(asctime)s [%(levelname)s] %(message)s"
		formatter = logging.Formatter(logFormat)
		ch.setFormatter(formatter)
		fh.setFormatter(formatter)
		
		# add the handlers to logger
		self.logger.addHandler(ch)
		self.logger.addHandler(fh)

		logging.info("========== Configuration ===========")
		logging.info("logFile:                 %s" % self.opts.logFile)
		logging.info("debug:                   %s" % self.opts.debug)
		logging.info("info:                    %s" % self.opts.info)
		logging.info("motif:                   %s" % self.opts.motif)
		logging.info("mod_pos:                 %s" % self.opts.mod_pos)
		logging.info("procs:                   %s" % self.opts.procs)
		logging.info("natProcs:                %s" % self.opts.natProcs)
		logging.info("chunkSize:               %s" % self.opts.chunkSize)
		logging.info("profile:                 %s" % self.opts.profile)
		logging.info("out:                     %s" % self.opts.out)
		logging.info("nativeCovThresh:         %s" % self.opts.nativeCovThresh)
		logging.info("wgaCovThresh:            %s" % self.opts.wgaCovThresh)
		logging.info("align:                   %s" % self.opts.align)
		logging.info("firstBasesToSkip:        %s" % self.opts.firstBasesToSkip)
		logging.info("lastBasesToSkip:         %s" % self.opts.lastBasesToSkip)
		logging.info("noCompare:               %s" % self.opts.noCompare)
		logging.info("nat_aligns_flat:         %s" % self.opts.nat_aligns_flat)
		logging.info("wga_aligns_flat:         %s" % self.opts.wga_aligns_flat)
		logging.info("extractAlignmentOnly:    %s" % self.opts.extractAlignmentOnly)
		logging.info("upstreamSkip:            %s" % self.opts.upstreamSkip)
		logging.info("downstreamSkip:          %s" % self.opts.downstreamSkip)
		logging.info("minSubreadLen:           %s" % self.opts.minSubreadLen)
		logging.info("minAcc:                  %s" % self.opts.minAcc)
		logging.info("minMapQV:                %s" % self.opts.minMapQV)
		logging.info("contig_name:             %s" % self.opts.contig_name)
		logging.info("SMsn:                    %s" % self.opts.SMsn)
		logging.info("SMp:                     %s" % self.opts.SMp)
		logging.info("wga_lib:                 %s" % self.opts.wga_lib)
		logging.info("nat_lib:                 %s" % self.opts.nat_lib)
		logging.info("leftAnchor:              %s" % self.opts.leftAnchor)
		logging.info("rightAnchor:             %s" % self.opts.rightAnchor)
		logging.info("write_vars:              %s" % self.opts.write_vars)
		logging.info("====================================")

	def get_reference_contigs( self, h5file ):
		"""
		Pull out the list of contigs in the h5 file.
		"""
		contigs   = map(lambda x: x.strip("/"), list(h5file["/RefGroup/Path"].value))
		for contig in contigs:
			logging.debug("  %s" % contig)
		return contigs

	def get_movie_name_ID_map( self, h5file ):
		"""
		Pull out the movie names and IDs from the h5 file and return a dict mapping them.
		"""
		movie_name_ID_map = dict(zip(h5file["/MovieInfo/Name"].value, h5file["/MovieInfo/ID"].value))
		for name, ID in movie_name_ID_map.iteritems():
			logging.debug("  %s : %s" % (name, ID))
		return movie_name_ID_map

	def extract_alignments_from_cmph5( self, cmph5, prefix ):
		"""
		Get movie names/IDs and contig names, then divide alignments in the aligned_reads.cmp.h5
		file into chunks that we load into memory. Once loaded as subread/alignment objects, we'll
		process them and write to the appropriate IPD csv file.
		"""
		h5file = h5py.File(cmph5)

		logging.debug("Getting all the reference contigs...")
		contigs = self.get_reference_contigs( h5file )

		logging.debug("Getting all the movie names and movie IDs...")
		movie_name_ID_map = self.get_movie_name_ID_map( h5file )
		
		# Too many chunks causes some I/O problems when reading from the cmp.h5
		nchunks         = min(self.opts.procs, 1)
		
		def write_chunks_of_alignments( chunk_id, chunk ):
			fn = multiprocessing.current_process().name
			f = open(fn, "w")
			tmp_flat_files.append(fn)
			n_IO_failures = 0
			for i, alignment in enumerate(chunk):
				if i%1000==0:
					logging.info("...chunk %s - alignment %s/%s (%.1f%%)" % (chunk_id, i, len(chunk), 100*float(i)/len(chunk)))
				try:
					if alignment.accuracy < self.opts.minAcc or alignment.MapQV < self.opts.minMapQV or len(alignment.alignmentArray()) < self.opts.minSubreadLen:
						continue
					movieID       = str(alignment.movieInfo[0])
					align_len     = str(len(alignment.alignmentArray()))
					fps           = str(alignment.movieInfo[2])
					refName       = str(alignment.referenceInfo[2])
					zmw           = str(alignment.HoleNumber)
					mol           = str(alignment.MoleculeID)
					if alignment.isForwardStrand:
						strand = 0
					else:
						strand = 1
					ref_bases  = alignment.reference()
					read_calls = alignment.transcript()
					ref_pos    = list(alignment.referencePositions())
					IPD        = list(alignment.IPD())
				except IOError:
					n_IO_failures += 1
					n_tolerated    = int(math.ceil(len(chunk) * 0.001))
					if n_IO_failures > n_tolerated:
						raise Exception("Too many IOErrors! Aborting!")
					logging.info("*** IOError in %s: %s so far. %s tolerated (0.1%% of all alignments) ***" % (fn, n_IO_failures, n_tolerated))
					continue

				if strand==0:
					f.write("%s %s %s %s %s %s %s %s %s %s %s\n" % (movieID, \
															 align_len, \
															 fps, \
															 refName, \
															 zmw, \
															 mol, \
															 strand, \
															 ref_bases, \
															 read_calls, \
															 ",".join(map(lambda x: str(x), ref_pos)), \
															 ",".join(map(lambda x: str(x), IPD))))
				elif strand==1:
					f.write("%s %s %s %s %s %s %s %s %s %s %s\n" % (movieID, \
															 align_len, \
															 fps, \
															 refName, \
															 zmw, \
															 mol, \
															 strand, \
															 ref_bases[::-1], \
															 read_calls[::-1], \
															 ",".join(map(lambda x: str(x), ref_pos[::-1])), \
															 ",".join(map(lambda x: str(x), IPD[::-1]))))
			f.close()

		def chunks(l, n):
			"""
			Yield successive n-sized chunks from l.
			"""
			for i in xrange(0, len(l), n):
				yield l[i:i+n]

		logging.info("Loading %s..." % cmph5)
		reader           = CmpH5Reader(cmph5)
		alignments_list  = [r for r in reader]
		logging.info("Done.")
		if (prefix == "nat_" and self.opts.nat_aligns_flat == None) or (prefix == "wga_" and self.opts.wga_aligns_flat == None):
			chunksize        = int(math.ceil(float( len(alignments_list)/nchunks )))
			alignment_chunks = list(chunks( alignments_list, chunksize ))
			workers          = []
			tmp_flat_files   = []
			for i,chunk in enumerate(alignment_chunks):
				worker_name = "%saligned_reads.txt.chunk%s" % (prefix, i)
				worker      = multiprocessing.Process(name=worker_name, target=write_chunks_of_alignments, args=(i, chunk,) )
				workers.append(worker)
				
			for worker in workers:
				worker.start()

			tmp_flat_files = []
			for worker in workers:
				worker.join()
				tmp_flat_files.append(worker.name)

			# catted_flat_file = "%saligned_reads.txt" % prefix
			
			catted_flat_file = cmph5.replace("cmp.h5", "flat.txt")
			catted_flat_file = cat_list_of_files( tmp_flat_files, catted_flat_file )
		elif prefix == "nat_":
			catted_flat_file = self.opts.nat_aligns_flat
		elif prefix == "wga_":
			catted_flat_file = self.opts.wga_aligns_flat

		return catted_flat_file, movie_name_ID_map

	def sort_flat_alignments_file( self, flat_file ):
		"""
		"""
		if os.path.exists("%s.sorted" % flat_file):
			pass
		else:
			sort_CMD  = "sort -t' ' -gk10 %s > %s.sorted" % (flat_file, flat_file)
			p         = subprocess.Popen(sort_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdOutErr = p.communicate()
			sts       = p.returncode
			if sts != 0:
				raise Exception("Failed system command: %s" % sort_CMD)
		return "%s.sorted" % flat_file

	def split_flat_file_by_nprocs( self, flat_file, prefix ):
		"""
		To facilitate multiprocessing, split up flat alignments file into <procs> separate files.
		"""
		wc_CMD    = "wc -l %s" % flat_file
		p         = subprocess.Popen(wc_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdOutErr = p.communicate()
		sts       = p.returncode
		if sts != 0:
			raise Exception("Failed system command: %s" % wc_CMD)
		nlines    = int(stdOutErr[0].split()[0])

		if prefix == "nat_":
			lines_per_file = int(math.ceil( float(nlines) / self.opts.natProcs ))
		else:
			lines_per_file = int(math.ceil( float(nlines) / self.opts.procs ))

		split_CMD      = "split -a 2 -l %s %s %s." % (lines_per_file, flat_file, flat_file)
		p              = subprocess.Popen(split_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdOutErr      = p.communicate()
		sts            = p.returncode
		if sts != 0:
			raise Exception("Failed system command: %s" % split_CMD)

		files = glob.glob("%s.*" % flat_file)
		files.sort()

		# Reunite any subreads from a molecule that have been split up into separate files
		mol_file_map   = {}
		all_split_mols = []
		for i,fn in enumerate(files):
			file_mols = set()
			for line in open(fn).xreadlines():
				mol_id = int(line.split(" ")[5])
				file_mols.add(mol_id)
			mol_file_map[fn] = file_mols

			if i>0:
				u = set.intersection( mol_file_map[files[i]], mol_file_map[files[i-1]] )
				all_split_mols += list(u)
		all_split_mols = set(all_split_mols)

		if prefix == "nat_":
			if len(files) != self.opts.natProcs:
				raise Exception("The number of split native alignment files doesn't match natProcs!")
		else:
			if len(files) != self.opts.procs:
				raise Exception("The number of split WGA alignment files doesn't match procs!")
		return files, all_split_mols

	def split_up_control_IPDs( self, control_ipds, tmp_flat_files ):
		"""
		Separate out relevant portions of the control_ipds dictionary. We are taking 
		advantage of the fact that the alignment flat files are sorted by aligned 
		reference position.
		"""
		def pull_last_ref_pos_from_alignments_file( fn, head_or_tail ):
			tail_CMD   = "%s -1 %s" % (head_or_tail, fn)
			p         = subprocess.Popen(tail_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdOutErr = p.communicate()
			sts       = p.returncode
			if sts != 0:
				raise Exception("Failed system command: %s" % tail_CMD)
			line = stdOutErr[0]
			if head_or_tail == "head":
				value = int(line.split(" ")[-2].split(",")[1])
			else:
				value = int(line.split(" ")[-2].split(",")[-1])
			return value

		local_control_ipds      = {}
		for chunk_id,alignments_flat_fn in enumerate(tmp_flat_files):
			first_ref_pos  = pull_last_ref_pos_from_alignments_file( alignments_flat_fn, "head" )
			last_ref_pos   = pull_last_ref_pos_from_alignments_file( alignments_flat_fn, "tail" )
			region_control = {0:{}, 1:{}}
			logging.debug("Split control IPD dicts  --  chunk %s: %sbp - %sbp" % (chunk_id, first_ref_pos, last_ref_pos+1))
			for strand in region_control.keys():
				for pos in range(first_ref_pos, last_ref_pos+1):
					try:
						region_control[strand][pos] = control_ipds[strand][pos]
					except KeyError:
						# In case we don't have WGA coverage at this position
						pass
			local_control_ipds[chunk_id]      = region_control
		return local_control_ipds

	def find_motif_sites(self):
		"""
		This will write two files for a given motif, modification position in the motif, 
		and reference fasta. One file will have the motif positions on the forward strand
		and one on the negative strand of the reference.
		"""
		caller_script    = os.path.join(os.path.dirname(__file__),'R/call_motifFinder.r')
		findMotif_script = os.path.join(os.path.dirname(__file__),'R/motifFinder.r')
		Rscript_CMD      = "Rscript %s %s %s %s %s" % (caller_script, findMotif_script, self.opts.motif, self.opts.mod_pos, self.ref)
		p                = subprocess.Popen(Rscript_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdOutErr        = p.communicate()
		sts              = p.returncode
		if sts != 0:
			for entry in stdOutErr:
				print entry
			raise Exception("Failed command: %s" % Rscript_CMD)
		self.sites_pos = stdOutErr[0].split("\n")[0].split(" ")[1][1:-1]
		self.sites_neg = stdOutErr[0].split("\n")[1].split(" ")[1][1:-1]

	def launch_parallel_molecule_loading( self, tmp_flat_files, prefix, movie_name_ID_map, control_ipds, split_mols ):
		"""
		"""
		if prefix == "nat_":
			logging.info("Separating out file-matched regions of the control IPD values dict...")
			local_control_ipds = self.split_up_control_IPDs( control_ipds, tmp_flat_files )
			logging.info("Done.")

		logging.debug("Creating tasks...")
		tasks   = multiprocessing.JoinableQueue()
		results = multiprocessing.Queue()
		logging.debug("Done.")

		if prefix == "nat_":
			num_consumers = self.opts.natProcs
		else:
			num_consumers = self.opts.procs

		logging.debug("Starting consumers...")
		consumers     = [ Consumer(tasks, results) for i in xrange(num_consumers) ]
		for w in consumers:
			w.start()
		logging.debug("Done.")

		# Enqueue jobs
		num_jobs = len(tmp_flat_files)
		for chunk_id,alignments_flat_fn in enumerate(tmp_flat_files):
			if prefix == "wga_":
				tasks.put(parse_mol_aligns.wga_molecules_processor( alignments_flat_fn,    \
																	chunk_id,              \
																	prefix,                \
																	self.opts.contig_name,  \
																	self.ref_size,         \
																	self.sites_pos,        \
																	self.sites_neg,        \
																	self.opts.leftAnchor,  \
																	self.opts.rightAnchor, \
																	self.opts.wga_lib,     \
																	split_mols ))
			else:
				logging.debug("Launching subprocess %s..." % chunk_id)
				tasks.put(parse_mol_aligns.native_molecules_processor( alignments_flat_fn,           \
																	   chunk_id,                     \
																	   prefix,                       \
																	   self.opts.contig_name,         \
																	   self.opts.nativeCovThresh,    \
																	   self.fastq,                   \
																	   self.ref,                     \
																	   self.opts.align,              \
																	   copy.copy(movie_name_ID_map), \
																	   self.opts.firstBasesToSkip,   \
																	   self.opts.lastBasesToSkip,    \
																	   self.opts.upstreamSkip,       \
																	   self.opts.downstreamSkip,     \
																	   local_control_ipds[chunk_id], \
																	   self.ref_size,                \
																	   self.sites_pos,               \
																	   self.sites_neg,               \
																	   self.opts.SMp,                \
																	   self.opts.leftAnchor,         \
																	   self.opts.rightAnchor,        \
																	   self.opts.nat_lib,            \
																	   split_mols,                   \
																	   self.opts.wgaCovThresh,       \
																	   self.opts.out))
				logging.debug("Done (%s)." % chunk_id)
		
		# Add a 'poison pill' for each consumer
		for i in xrange(num_consumers):
			tasks.put(None)
		tasks.join()
		
		# Start printing results
		parallel_results = []
		while num_jobs:
			result = results.get()
			parallel_results.append(result)
			num_jobs -= 1

		return parallel_results

	def run( self ):
		"""
		Execute the pipeline.
		"""
		self.find_motif_sites()

		##############
		# WGA
		##############
		prefix = "wga_"
		logging.info("Extracting WGA alignments from %s..." % self.wga_cmph5)
		wga_flat_file, wga_movie_name_ID_map = self.extract_alignments_from_cmph5( self.wga_cmph5, prefix )
		logging.info("Done.")

		logging.info("Sorting %s..." % wga_flat_file)
		wga_flat_file = self.sort_flat_alignments_file( wga_flat_file )
		logging.info("Done.")

		if not self.opts.extractAlignmentOnly:
			logging.info("Splitting %s into %s chunks for multiprocessing..." % (wga_flat_file, self.opts.procs))
			wga_tmp_flat_files, split_mols = self.split_flat_file_by_nprocs( wga_flat_file, prefix )
			logging.info("Done.")

			control_ipds    = None
			chunk_ipdArrays = self.launch_parallel_molecule_loading( wga_tmp_flat_files, prefix, wga_movie_name_ID_map, control_ipds, split_mols )

			logging.info("Combining the %s separate ipdArray dictionaries..." % len(chunk_ipdArrays))
			control_ipds    = parse_mol_aligns.combine_chunk_ipdArrays( chunk_ipdArrays, self.ref_size )
			logging.info("Done.")

			logging.info("Cleaning up chunked WGA alignment files...")
			for fn in wga_tmp_flat_files:
				os.remove(fn)
			logging.info("Done")
			
		##############
		# Native
		##############
		prefix = "nat_"
		logging.info("Extracting native alignments from %s..." % self.native_cmph5)
		native_flat_file, native_movie_name_ID_map = self.extract_alignments_from_cmph5( self.native_cmph5, prefix )
		logging.info("Done.")
		
		logging.info("Sorting %s..." % native_flat_file)
		native_flat_file = self.sort_flat_alignments_file( native_flat_file )
		logging.info("Done.")

		if not self.opts.extractAlignmentOnly:
			logging.info("Splitting %s into %s chunks for multiprocessing..." % (native_flat_file, self.opts.natProcs))
			native_tmp_flat_files, split_mols = self.split_flat_file_by_nprocs( native_flat_file, prefix )
			logging.info("Done.")

			if self.opts.align:
				# Create necessary fasta index for BWA aligner and samtools
				logging.info("Indexing %s..." % self.ref)
				bwa_idx_CMD      = "bwa index %s" % (self.ref)
				samtools_idx_CMD = "samtools faidx %s" % self.ref
				run_command( bwa_idx_CMD )
				run_command( samtools_idx_CMD )
				logging.info("Done.")

			parallel_output_fns = self.launch_parallel_molecule_loading( native_tmp_flat_files, prefix, native_movie_name_ID_map, control_ipds, split_mols )

			logging.info("Cleaning up chunked native alignment files...")
			for fn in native_tmp_flat_files:
				os.remove(fn)
			logging.info("Done.")

			logging.info("Combining chunked test output files...")
			out_files_to_cat = [fn for fn in parallel_output_fns if os.path.exists(fn)]
			head             = "strand\tpos\tscore\tmol\tnat\twga\tN_nat\tN_wga\tsubread_len\n"
			self.opts.out    = cat_list_of_files( out_files_to_cat, self.opts.out, header=head )
			logging.info("Done.")

			if self.opts.align and self.opts.write_vars:
				logging.info("Combining chunked molecule-specific variant calls...")
				vars_files_to_cat = glob.glob("vars_*.tmp")
				head = "mol\tvar_pos\n"
				self.opts.write_vars = cat_list_of_files(vars_files_to_cat , self.opts.write_vars, header=head )
				logging.info("Done.")
			elif self.opts.align:
				vars_files_to_del = glob.glob("vars_*.tmp")
				for fn in vars_files_to_del:
					os.remove(fn)

			logging.info("Finalizing cleanup...")
			for fn in glob.glob("%s.*" % self.ref):
				os.remove(fn)
			for fn in glob.glob("*flat.txt.sorted"):
				os.remove(fn)
			logging.info("Done.")

def main():
	app = Smalr_runner()
	if app.opts.profile:
		import cProfile
		cProfile.run( 'app.run()')
	else:
		sys.exit( app.run() )