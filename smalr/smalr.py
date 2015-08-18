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
	# for fn in in_fns:
	# 	try:
	# 		os.remove(fn)
	# 	except OSError:
	# 		pass
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

	def get_movie_name_ID_map( self, h5file ):
		"""
		Pull out the movie names and IDs from the h5 file and return a dict mapping them.
		"""
		movie_name_ID_map = dict(zip(h5file["/MovieInfo/Name"].value, h5file["/MovieInfo/ID"].value))
		for name, ID in movie_name_ID_map.iteritems():
			logging.debug("  %s : %s" % (name, ID))
		return movie_name_ID_map

	def split_up_control_IPDs( self, control_ipds, cmph5_file, idx_chunks ):
		"""
		Separate out relevant portions of the control_ipds dictionary. We are taking 
		advantage of the fact that the alignment flat files are sorted by aligned 
		reference position.
		"""

		reader                  = CmpH5Reader(cmph5_file)
		local_control_ipds      = {}
		for chunk_id,idx_chunk in enumerate(idx_chunks):
			idx_mins = [min(reader[idx].tStart, reader[idx].tEnd) for idx in idx_chunk]
			idx_maxs = [max(reader[idx].tStart, reader[idx].tEnd) for idx in idx_chunk]
			first_ref_pos  = min(idx_mins)
			last_ref_pos   = max(idx_maxs)
			# first_ref_pos  = pull_last_ref_pos_from_alignments_file( alignments_flat_fn, "head" )
			# last_ref_pos   = pull_last_ref_pos_from_alignments_file( alignments_flat_fn, "tail" )
			region_control = {0:{}, 1:{}}
			logging.debug("Split control IPD dicts  --  chunk %s: %sbp - %sbp" % (chunk_id, first_ref_pos, last_ref_pos+1))
			for strand in region_control.keys():
				for pos in range(first_ref_pos, last_ref_pos+1):
					try:
						region_control[strand][pos] = control_ipds[strand][pos]
					except KeyError:
						# In case we don't have WGA coverage at this position
						pass
			local_control_ipds[chunk_id] = region_control
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

		if not os.path.exists(contig_fasta_fn) or os.path.getsize(contig_fasta_fn)==0:
			raise Exception("Couldn't write the contig-specific fasta file!")

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

	def track_split_molecule_alignments( self, idx_chunks, cmph5_file ):
		reader     = CmpH5Reader(cmph5_file)
		chunk_mols = {}
		for i,idx_chunk in enumerate(idx_chunks):
			chunk_mols[i] = set()
			for alignment in reader[idx_chunk]:
				chunk_mols[i].add(alignment.MoleculeID)
		reader.close()

		split_mols = set()
		i = 1
		for idx_chunk in idx_chunks[1:]:
			j = i-1
			split = chunk_mols[i] & chunk_mols[j]
			split_mols = split_mols | split
			i += 1
		return split_mols

	def launch_parallel_molecule_loading( self, cmph5_file, prefix, movie_name_ID_map, control_ipds ):
		logging.debug("Creating tasks...")
		tasks   = multiprocessing.JoinableQueue()
		results = multiprocessing.Queue()
		logging.debug("Done.")

		logging.debug("Starting consumers...")
		num_jobs      = self.Config.opts.procs
		consumers     = [ Consumer(tasks, results, self.Config.opts.contig_id) for i in xrange(num_jobs) ]
		for w in consumers:
			w.start()
		logging.debug("Done.")

		def chunks( l, n ):
			"""
			Yield successive n-sized chunks from l.
			"""
			for i in xrange(0, len(l), n):
				yield l[i:i+n]

		# Enqueue jobs
		logging.info("Partitioning %s into %s chunks for analysis..." % (cmph5_file, num_jobs))
		reader         = CmpH5Reader(cmph5_file)
		alnIDs         = [r.AlnID for r in reader if r.referenceInfo[2]==self.Config.opts.contig_id]
		if len(alnIDs) <= num_jobs:
			num_jobs = 1
		reader.close()

		# for chunk_id,alignments_flat_fn in enumerate(tmp_flat_files):
		chunksize    = int(math.ceil(float( len(alnIDs)/num_jobs )))
		idx_chunks   = list(chunks( (np.array(alnIDs)-1), chunksize ))

		if len(idx_chunks[-1])==1:
			idx_chunks = idx_chunks[:-1]

		if prefix == "nat_":
			logging.info("%s - Separating out file-matched regions of the control IPD values dict..." % self.Config.opts.contig_id)
			local_control_ipds = self.split_up_control_IPDs( control_ipds, cmph5_file, idx_chunks )
			logging.debug("%s - Done." % self.Config.opts.contig_id)

		# In splitting alignment indexes among processes, some molecules will have 
		# alignments in going to different processes. Track these.
		split_mols = self.track_split_molecule_alignments( idx_chunks, cmph5_file )

		for chunk_id in range(num_jobs):
			idx = idx_chunks[chunk_id]
			if prefix == "wga_":
				tasks.put(parse_mol_aligns.wga_molecules_processor( cmph5_file,                   \
																	chunk_id,                     \
																	prefix,                       \
																	self.Config.opts.contig_id,   \
																	self.ref_size,                \
																	self.sites_pos,               \
																	self.sites_neg,               \
																	self.Config.opts,             \
																	idx,                          \
																	split_mols ))
			else:
				logging.debug("Launching subprocess %s..." % chunk_id)
				tasks.put(parse_mol_aligns.native_molecules_processor( cmph5_file,                          \
																	   chunk_id,                            \
																	   prefix,                              \
																	   self.Config.opts,                    \
																	   self.Config.fastq,                   \
																	   self.Config.ref,                     \
																	   copy.copy(movie_name_ID_map),        \
																	   local_control_ipds[chunk_id],        \
																	   self.ref_size,                       \
																	   self.sites_pos,                      \
																	   self.sites_neg,                      \
																	   idx,                                 \
																	   split_mols))
				logging.debug("Done (%s)." % chunk_id)
		
		# Add a 'poison pill' for each consumer
		for i in xrange(num_jobs):
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
		wga_movie_name_ID_map = self.get_movie_name_ID_map( h5file = h5py.File(self.Config.wga_cmph5, 'r') )

		control_ipds    = None
		chunk_ipdArrays = self.launch_parallel_molecule_loading( self.Config.wga_cmph5, prefix, wga_movie_name_ID_map, control_ipds )

		logging.info("%s - Combining the %s separate ipdArray dictionaries..." % (self.Config.opts.contig_id, len(chunk_ipdArrays)))
		control_ipds    = parse_mol_aligns.combine_chunk_ipdArrays( chunk_ipdArrays, self.ref_size )
		logging.debug("%s - Done." % self.Config.opts.contig_id)
			
		##############
		# Native
		##############
		prefix = "nat_"
		if self.Config.opts.align:
			# Create necessary fasta index for BWA aligner and samtools
			logging.info("%s - Indexing %s..." % (self.Config.opts.contig_id, self.Config.ref))
			bwa_idx_CMD      = "bwa index %s" % (self.Config.ref)
			samtools_idx_CMD = "samtools faidx %s" % self.Config.ref
			run_command( bwa_idx_CMD )
			run_command( samtools_idx_CMD )
			logging.debug("%s - Done." % self.Config.opts.contig_id)

		native_movie_name_ID_map = self.get_movie_name_ID_map( h5file = h5py.File(self.Config.native_cmph5, 'r') )

		parallel_output_fns = self.launch_parallel_molecule_loading( self.Config.native_cmph5, prefix, native_movie_name_ID_map, control_ipds )

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
		logging.debug("%s - Done." % self.Config.opts.contig_id)

		if self.Config.opts.out != None:
			logging.info("Sorting output by strand/position...")
			sort_CMD = "sort -t$'\\t' -nk2 %s > sorting.tmp" % self.Config.opts.out
			p         = subprocess.Popen(sort_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdOutErr = p.communicate()
			sts       = p.returncode
			if sts != 0:
				raise Exception("Failed sort command: %s" % sort_CMD)
			os.rename("sorting.tmp", self.Config.opts.out)
			logging.info("Done.")

def main():
	app = SmalrRunner()
	sys.exit( app.run() )
