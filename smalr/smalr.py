import os,sys,time
import parse_mol_aligns
import subprocess
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
import copy
import operator
import SmalrConfig

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
		return

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

class Consumer(multiprocessing.Process):
	def __init__(self, task_queue, result_queue, contig_id):
		multiprocessing.Process.__init__(self)
		self.task_queue   = task_queue
		self.result_queue = result_queue
		self.contig_id    = contig_id

	def run(self):
		proc_name = self.name
		while True:
			next_task = self.task_queue.get()
			if next_task is None:
				# Poison pill means shutdown
				logging.debug("%s - %s: Exiting" % (self.contig_id, proc_name))
				self.task_queue.task_done()
				break
			logging.debug("%s - %s: Starting" % (self.contig_id, proc_name))
			answer = next_task()
			self.task_queue.task_done()
			self.result_queue.put(answer)
		return

class SmalrRunner():
	def __init__ ( self, i, contig_info, abs_input_file, Config ):
		"""
		Parse the options and arguments, then instantiate the logger. 
		"""
		self.i                       = i
		self.Config                  = Config
		self.Config.opts.contig_name = contig_info[0]
		self.Config.opts.contig_id   = contig_info[1]
		self.__initLog( )
		
		logging.info("%s - contig_id:               %s" % (self.Config.opts.contig_id, self.Config.opts.contig_id))
		logging.info("%s - contig_name:             %s" % (self.Config.opts.contig_id, self.Config.opts.contig_name))

		for i,entry in enumerate(fasta_iter(self.Config.ref)):
			self.ref_name = entry[0]
			if self.ref_name == self.Config.opts.contig_name:
				self.ref_size = len(entry[1])
			else:
				pass

	def __initLog( self ):
		"""Sets up logging based on command line arguments. Allows for three levels of logging:
		logging.error( ): always emitted
		logging.info( ) : emitted with --info or --debug
		logging.debug( ): only with --debug"""

		if os.path.exists(self.Config.opts.logFile):
			os.remove(self.Config.opts.logFile)

		logLevel = logging.DEBUG if self.Config.opts.debug \
					else logging.INFO if self.Config.opts.info \
					else logging.ERROR

		self.logger = logging.getLogger()
		self.logger.setLevel(logLevel)

		# create file handler which logs even debug messages
		fh = logging.FileHandler(self.Config.opts.logFile)
		fh.setLevel(logLevel)
		
		# create console handler with a higher log level
		ch = logging.StreamHandler()
		ch.setLevel(logLevel)
		
		# create formatter and add it to the handlers
		logFormat = "%(asctime)s [%(levelname)s] %(message)s"
		formatter = logging.Formatter(logFormat, datefmt='%H:%M:%S')
		ch.setFormatter(formatter)
		fh.setFormatter(formatter)
		
		# add the handlers to logger
		if self.i == 0:
			self.logger.addHandler(ch)
		self.logger.addHandler(fh)
		logging.info("")
		logging.info("====================================")
		logging.info("Analyzing contig %s (%s)" % (self.Config.opts.contig_id, self.Config.opts.contig_name))
		logging.info("====================================")

	def get_reference_contigs( self, h5file ):
		"""
		Pull out the list of contigs in the h5 file.
		"""
		contigs   = map(lambda x: x.strip("/"), list(h5file["/RefGroup/Path"].value))
		for contig in contigs:
			if contig == self.Config.opts.contig_id:
				logging.info("%s -   %s  <-- Analyzing this (%s)" % (self.Config.opts.contig_id, contig, self.Config.opts.contig_name))
			else:
				logging.info("%s -   %s" % (self.Config.opts.contig_id, contig))
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
		h5file = h5py.File(cmph5, 'r')
		logging.debug("Getting all the reference contigs...")
		contigs = self.get_reference_contigs( h5file )
		logging.debug("Getting all the movie names and movie IDs...")
		movie_name_ID_map = self.get_movie_name_ID_map( h5file )
		h5file.close()
		
		# Too many chunks causes some I/O problems when reading from the cmp.h5
		nchunks         = min(self.Config.opts.procs, 1)
		
		def write_chunks_of_alignments( chunk_id, chunk ):
			fn = multiprocessing.current_process().name
			f = open(fn, "w")
			tmp_flat_files.append(fn)
			n_IO_failures = 0
			for i, alignment in enumerate(chunk):
				if i%1000==0:
					logging.info("%s - ...chunk %s - alignment %s/%s (%.1f%%)" % (self.Config.opts.contig_id,chunk_id, i, len(chunk), 100*float(i)/len(chunk)))
				try:
					if alignment.accuracy < self.Config.opts.minAcc or alignment.MapQV < self.Config.opts.minMapQV or len(alignment.alignmentArray()) < self.Config.opts.minSubreadLen:
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
					logging.info("%s - *** IOError in %s: %s so far. %s tolerated (0.1%% of all alignments) ***" % (self.Config.opts.contig_id, fn, n_IO_failures, n_tolerated))
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

		logging.info("%s - Loading %s..." % (self.Config.opts.contig_id, cmph5))
		reader           = CmpH5Reader(cmph5)
		alignments_list  = [r for r in reader]
		logging.debug("%s - Done." % self.Config.opts.contig_id)
		if (prefix == "nat_" and self.Config.opts.nat_aligns_flat == None) or (prefix == "wga_" and self.Config.opts.wga_aligns_flat == None):
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
			catted_flat_file = self.Config.opts.nat_aligns_flat
		elif prefix == "wga_":
			catted_flat_file = self.Config.opts.wga_aligns_flat

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
			lines_per_file = int(math.ceil( float(nlines) / self.Config.opts.natProcs ))
		else:
			lines_per_file = int(math.ceil( float(nlines) / self.Config.opts.procs ))

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
			if len(files) != self.Config.opts.natProcs:
				raise Exception("The number of split native alignment files doesn't match natProcs!")
		else:
			if len(files) != self.Config.opts.procs:
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
		f_iter = fasta_iter(self.Config.ref)
		
		contig_fasta_fn = "%s.fasta" % self.Config.opts.contig_id
		f               = open(contig_fasta_fn, "w")
		for name,seq in f_iter:
			if name == self.Config.opts.contig_name:
				f.write(">%s\n" % name)
				f.write("%s" % seq)
		f.close()

		caller_script    = os.path.join(os.path.dirname(__file__),'R/call_motifFinder.r')
		findMotif_script = os.path.join(os.path.dirname(__file__),'R/motifFinder.r')
		Rscript_CMD      = "Rscript %s %s %s %s %s" % (caller_script, findMotif_script, self.Config.opts.motif, self.Config.opts.mod_pos, contig_fasta_fn)
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
			logging.info("%s - Separating out file-matched regions of the control IPD values dict..." % self.Config.opts.contig_id)
			local_control_ipds = self.split_up_control_IPDs( control_ipds, tmp_flat_files )
			logging.debug("%s - Done." % self.Config.opts.contig_id)

		logging.debug("Creating tasks...")
		tasks   = multiprocessing.JoinableQueue()
		results = multiprocessing.Queue()
		logging.debug("Done.")

		if prefix == "nat_":
			num_consumers = self.Config.opts.natProcs
		else:
			num_consumers = self.Config.opts.procs

		logging.debug("Starting consumers...")
		consumers     = [ Consumer(tasks, results, self.Config.opts.contig_id) for i in xrange(num_consumers) ]
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
																	self.Config.opts.contig_id,  \
																	self.ref_size,         \
																	self.sites_pos,        \
																	self.sites_neg,        \
																	self.Config.opts.leftAnchor,  \
																	self.Config.opts.rightAnchor, \
																	self.Config.opts.wga_lib,     \
																	split_mols ))
			else:
				logging.debug("Launching subprocess %s..." % chunk_id)
				tasks.put(parse_mol_aligns.native_molecules_processor( alignments_flat_fn,           \
																	   chunk_id,                     \
																	   prefix,                       \
																	   self.Config.opts.contig_id,         \
																	   self.Config.opts.nativeCovThresh,    \
																	   self.Config.fastq,                   \
																	   self.Config.ref,                     \
																	   self.Config.opts.align,              \
																	   copy.copy(movie_name_ID_map), \
																	   self.Config.opts.firstBasesToSkip,   \
																	   self.Config.opts.lastBasesToSkip,    \
																	   self.Config.opts.upstreamSkip,       \
																	   self.Config.opts.downstreamSkip,     \
																	   local_control_ipds[chunk_id], \
																	   self.ref_size,                \
																	   self.sites_pos,               \
																	   self.sites_neg,               \
																	   self.Config.opts.SMp,                \
																	   self.Config.opts.leftAnchor,         \
																	   self.Config.opts.rightAnchor,        \
																	   self.Config.opts.nat_lib,            \
																	   split_mols,                   \
																	   self.Config.opts.wgaCovThresh,       \
																	   self.Config.opts.out))
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
		logging.info("%s - Extracting WGA alignments from %s..." % (self.Config.opts.contig_id, self.Config.wga_cmph5))
		wga_flat_file, wga_movie_name_ID_map = self.extract_alignments_from_cmph5( self.Config.wga_cmph5, prefix )
		logging.debug("%s - Done." % self.Config.opts.contig_id)

		logging.info("%s - Sorting %s..." % (self.Config.opts.contig_id, wga_flat_file))
		wga_flat_file = self.sort_flat_alignments_file( wga_flat_file )
		logging.debug("%s - Done." % self.Config.opts.contig_id)

		if not self.Config.opts.extractAlignmentOnly:
			logging.info("%s - Splitting %s into %s chunks for multiprocessing..." % (self.Config.opts.contig_id, wga_flat_file, self.Config.opts.procs))
			wga_tmp_flat_files, split_mols = self.split_flat_file_by_nprocs( wga_flat_file, prefix )
			logging.debug("%s - Done." % self.Config.opts.contig_id)

			control_ipds    = None
			chunk_ipdArrays = self.launch_parallel_molecule_loading( wga_tmp_flat_files, prefix, wga_movie_name_ID_map, control_ipds, split_mols )

			logging.info("%s - Combining the %s separate ipdArray dictionaries..." % (self.Config.opts.contig_id, len(chunk_ipdArrays)))
			control_ipds    = parse_mol_aligns.combine_chunk_ipdArrays( chunk_ipdArrays, self.ref_size )
			logging.debug("%s - Done." % self.Config.opts.contig_id)

			logging.info("%s - Cleaning up chunked WGA alignment files..." % self.Config.opts.contig_id)
			for fn in wga_tmp_flat_files:
				os.remove(fn)
			logging.debug("%s - Done" % self.Config.opts.contig_id)
			
		##############
		# Native
		##############
		prefix = "nat_"
		logging.info("%s - Extracting native alignments from %s..." % (self.Config.opts.contig_id, self.Config.native_cmph5))
		native_flat_file, native_movie_name_ID_map = self.extract_alignments_from_cmph5( self.Config.native_cmph5, prefix )
		logging.debug("%s - Done." % self.Config.opts.contig_id)
		
		logging.info("%s - Sorting %s..." % (self.Config.opts.contig_id, native_flat_file))
		native_flat_file = self.sort_flat_alignments_file( native_flat_file )
		logging.debug("%s - Done." % self.Config.opts.contig_id)

		if not self.Config.opts.extractAlignmentOnly:
			logging.info("%s - Splitting %s into %s chunks for multiprocessing..." % (self.Config.opts.contig_id, native_flat_file, self.Config.opts.natProcs))
			native_tmp_flat_files, split_mols = self.split_flat_file_by_nprocs( native_flat_file, prefix )
			logging.debug("%s - Done." % self.Config.opts.contig_id)

			if self.Config.opts.align:
				# Create necessary fasta index for BWA aligner and samtools
				logging.info("%s - Indexing %s..." % (self.Config.opts.contig_id, self.Config.ref))
				bwa_idx_CMD      = "bwa index %s" % (self.Config.ref)
				samtools_idx_CMD = "samtools faidx %s" % self.Config.ref
				run_command( bwa_idx_CMD )
				run_command( samtools_idx_CMD )
				logging.debug("%s - Done." % self.Config.opts.contig_id)

			parallel_output_fns = self.launch_parallel_molecule_loading( native_tmp_flat_files, prefix, native_movie_name_ID_map, control_ipds, split_mols )

			logging.info("%s - Cleaning up chunked native alignment files..." % self.Config.opts.contig_id)
			for fn in native_tmp_flat_files:
				os.remove(fn)
			logging.debug("%s - Done." % self.Config.opts.contig_id)

			logging.info("%s - Combining chunked test output files..." % self.Config.opts.contig_id)
			out_files_to_cat = [fn for fn in parallel_output_fns if os.path.exists(fn)]
			head             = "strand\tpos\tscore\tmol\tnat\twga\tN_nat\tN_wga\tsubread_len\n"
			self.Config.opts.out    = cat_list_of_files( out_files_to_cat, self.Config.opts.out, header=head )
			logging.debug("%s - Done." % self.Config.opts.contig_id)

			if self.Config.opts.align and self.Config.opts.write_vars:
				logging.info("%s - Combining chunked molecule-specific variant calls..." % self.Config.opts.contig_id)
				vars_files_to_cat = glob.glob("vars_*.tmp")
				head = "mol\tvar_pos\n"
				self.Config.opts.write_vars = cat_list_of_files(vars_files_to_cat , self.Config.opts.write_vars, header=head )
				logging.debug("%s - Done." % self.Config.opts.contig_id)
			elif self.Config.opts.align:
				vars_files_to_del = glob.glob("vars_*.tmp")
				for fn in vars_files_to_del:
					os.remove(fn)

			logging.info("%s - Finalizing cleanup..." % self.Config.opts.contig_id)
			for fn in glob.glob("%s.*" % self.Config.ref):
				os.remove(fn)
			for fn in glob.glob("*flat.txt.sorted"):
				os.remove(fn)
			logging.debug("%s - Done." % self.Config.opts.contig_id)

def main():
	app = SmalrRunner()
	sys.exit( app.run() )
