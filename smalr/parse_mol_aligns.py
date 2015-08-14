import sys,os
import shutil
import glob
import numpy as np
import math
from collections import defaultdict,Counter
from pbcore.io.align.CmpH5IO import *
import logging
import subprocess
import random
import copy
from operator import itemgetter
import CCS_aligner

class ipd_entry:
	def __init__( self, tup ):
		"""
		"""
		self.base    = tup[0]
		self.pos     = tup[1]
		self.call    = tup[2]
		self.ipd     = tup[3]
		self.strand  = np.abs(tup[4]-1)
		self.subread = tup[5]

class molecule:
	def __init__( self, alignments, prefix, leftAnchor, rightAnchor, contig_id, sites_pos, sites_neg ):
		"""
		"""
		self.leftAnchor  = leftAnchor
		self.rightAnchor = rightAnchor
		self.contig      = contig_id
		self.mapped      = False
		self.in_fastq    = False
		self.to_del      = False
		self.load_entries( alignments, sites_pos, sites_neg )
		
		# self.subread_normalize(  )

	def load_entries( self, alignments, sites_pos, sites_neg ):
		self.zmw_id        = "%s_%s" % (alignments[0].HoleNumber, alignments[0].MovieID)
		self.entries       = {}
		self.subread_means = {}
		for subread_id,alignment in enumerate(alignments):
			fps             = alignment.movieInfo[2]
			self.refName    = alignment.referenceInfo[3]
			self.mol_id     = alignment.MoleculeID
			if alignment.isForwardStrand:
				self.strand = 0
			else:
				self.strand = 1
			ref_bases  = alignment.reference()
			read_bases = alignment.read()
			read_calls = alignment.transcript()
			ref_pos    = list(alignment.referencePositions())
			IPD        = list(alignment.IPD())

			error_mk = []
			for read_call in read_calls:
				# Go through all entries and flag which positions are MM/indels
				if read_call != "M":
					# Mismatch or indel at this position!
					error_mk.append(1)
				else:
					error_mk.append(0)

			# Get the indices of all the non-matches
			error_idx = [i for (i,val) in enumerate(error_mk) if val == 1]
			for error_id in error_idx:
				try:
					for j in range(self.leftAnchor):
						error_mk[error_id - (j+1)] = 1
					for j in range(self.rightAnchor):
						error_mk[error_id + (j+1)] = 1
				except IndexError:
					pass
			error_mk = np.array(error_mk)

			ipds       = np.array(IPD) / fps
			strands    = np.array([self.strand] * len(read_calls))
			subread    = np.array([subread_id]  * len(read_calls))

			ref_bases  = np.array(list(ref_bases))
			read_bases = np.array(list(read_bases))
			ref_pos    = np.array(ref_pos)
			read_calls = np.array(list(read_calls))

			# Mark the error positions, but leave them in the sequence so
			# we can pull out intact motifs from contiguous correct bases
			ref_bases[error_mk==1]  = "*"
			read_bases[error_mk==1] = "*"
			read_calls[error_mk==1] = "*"
			ipds[error_mk==1]       = -9

			# Calculate subread mean IPD values
			if len(error_mk==0) < 2:
				self.subread_means[subread_id] = 0.1
			else:
				self.subread_means[subread_id] = np.log(ipds[error_mk==0] + 0.001).mean()

			if self.strand==1:
				keep_pos = list(set(ref_pos) & set(sites_pos))
			else:
				keep_pos = list(set(ref_pos) & set(sites_neg))

			motif_mk = []
			for pos in ref_pos[error_mk==0]:
				if pos in keep_pos:
					motif_mk.append(1)
				else:
					motif_mk.append(0)
			motif_mk = np.array(motif_mk)

			# Normalize these IPD entries and attach to the molecule object
			normed_ipds = np.exp(np.log(ipds[error_mk==0][motif_mk==1]+0.001) - self.subread_means[subread_id])
			for tup in zip(ref_bases[error_mk==0][motif_mk==1],   \
						   ref_pos[error_mk==0][motif_mk==1],     \
						   read_calls[error_mk==0][motif_mk==1],  \
						   normed_ipds,                           \
						   strands[error_mk==0][motif_mk==1],     \
						   subread[error_mk==0][motif_mk==1]):
				entry = ipd_entry(tup)
				self.entries[ (entry.pos, entry.subread) ] = entry

def read_in_motif_sites( motifSites_fn ):
	"""
	"""
	sites = []
	for line in open(motifSites_fn).xreadlines():
		sites.append(int(float(line.strip())))
	return sites

def generate_molecule_ZMW_map( mols, chunk_id ):
	"""
	Map zmw_id <--> mol_id. This will also remove any double-loaded molecules.
	"""
	zmw_mol_map = {}
	zmw_mol_counts   = Counter()
	for mol in mols.values():
		zmw_mol_map[mol.mol_id] = mol.zmw_id
		zmw_mol_counts[mol.zmw_id] += 1

	mol_ids              = zmw_mol_map.keys()
	multi_loaded_mol_ids = [mol_id for mol_id in mol_ids if zmw_mol_counts[zmw_mol_map[mol_id]] > 1]
	mols            = dict([ (mol.mol_id, mol) for mol in mols.values() if mol.mol_id not in multi_loaded_mol_ids])
	try:
		logging.debug("Process %s: removed %s double-loaded molecules (%.2f%%)." % (chunk_id,                  \
																			  		len(multi_loaded_mol_ids), \
																			  		100*len(multi_loaded_mol_ids)/float(len(mols.values()))))
	except ZeroDivisionError:
		logging.debug("Process %s: no molecules available for ZMW/moleculeID mapping!" % chunk_id)
	return mols, zmw_mol_map

def remove_split_up_molecules( mols, split_mols ):
	"""
	A few molecules have their alignments split between processes, so split_mols
	keeps track of these and now we'll exclude them from further analysis. This
	is a stop-gap measure until we can figure out a quick way of splitting the 
	original alignments file more cleverly.
	"""
	all_mols = set([mol.mol_id for mol in mols.values()])
	to_del   = list(set.intersection(split_mols, all_mols))
	for mol in to_del:
		del mols[mol]
	return mols

def remove_nonmotif_entries( ipdArrays, sites_pos, sites_neg, left=10, right=10 ):
	"""
	Go through the ipdArrays and only keep those that correspond to a motif site.
	"""
	ipdArrays_to_keep  = {0:{}, 1:{}}
	motif_pos_set_plus = set(sites_pos)
	ipds_pos_set_plus  = set(ipdArrays[0].keys())
	overlap_plus       = list(motif_pos_set_plus & ipds_pos_set_plus)
	if len(overlap_plus)>0:
		for p in overlap_plus:		
			ipdArrays_to_keep[0][p] = ipdArrays[0][p]

	motif_pos_set_minus = set(sites_neg)
	ipds_pos_set_minus  = set(ipdArrays[1].keys())
	overlap_minus       = list(motif_pos_set_minus & ipds_pos_set_minus)
	if len(overlap_minus)>0:
		for p in overlap_minus:		
			ipdArrays_to_keep[1][p] = ipdArrays[1][p]

	return ipdArrays_to_keep

def combine_chunk_ipdArrays( chunk_ipdArrays, ref_size ):
	"""
	"""
	combined_ipdArray = {0:{}, 1:{}}
	for ipdArray in chunk_ipdArrays:
		for strand in ipdArray.keys():
			for pos in ipdArray[strand].keys():
				try:
					combined_ipdArray[strand][pos] = np.concatenate( (ipdArray[strand][pos], combined_ipdArray[strand][pos]) )
				except KeyError:
					combined_ipdArray[strand][pos] = ipdArray[strand][pos]
	return combined_ipdArray

class wga_molecules_processor:
	def __init__( self,               \
				  cmph5,              \
				  chunk_id,           \
				  prefix,             \
				  contig_id,          \
				  ref_size,           \
				  sites_pos,          \
				  sites_neg,          \
				  opts,               \
				  idx,                \
				  split_mols ):
		"""
		"""
		self.cmph5              = cmph5
		self.chunk_id           = chunk_id
		self.prefix             = prefix
		self.contig_id          = contig_id
		self.ref_size           = ref_size
		self.sites_pos_fn       = sites_pos
		self.sites_neg_fn       = sites_neg
		self.leftAnchor         = opts.leftAnchor
		self.rightAnchor        = opts.rightAnchor
		self.wga_lib            = opts.wga_lib
		self.contig_id          = opts.contig_id
		self.opts               = opts
		self.idx                = idx
		self.split_mols         = split_mols

		if self.sites_pos_fn != None and self.sites_neg_fn != None:
			logging.debug("Process %s: reading motif sites..." % self.chunk_id)
			self.sites_pos = read_in_motif_sites( self.sites_pos_fn  )
			self.sites_neg = read_in_motif_sites( self.sites_neg_fn  )

	def __call__( self ):
		reader = CmpH5Reader(self.cmph5)

		mol_alignments = defaultdict(list)
		for alignment in reader[self.idx]:
			if alignment.accuracy >= self.opts.minAcc and alignment.MapQV >= self.opts.minMapQV and len(alignment.alignmentArray()) >= self.opts.minSubreadLen:
				mol_id = alignment.MoleculeID
				mol_alignments[mol_id].append(alignment)

		self.mols = {}
		i         = 0
		incr      = int(max(4,math.floor(float(len(mol_alignments.keys())))/4))
		for i,(mol_id, alignments) in enumerate(mol_alignments.iteritems()):
			if i%incr==0:
				logging.info("...chunk %s - %s/%s (%.1f%%) alignments processed..." % (self.chunk_id, i, len(mol_alignments.keys()), 100*float(i)/len(mol_alignments.keys())))
			mol_id                = alignments[0].MoleculeID
			mol                   = molecule( alignments, self.prefix, self.leftAnchor, self.rightAnchor, self.contig_id, self.sites_pos, self.sites_neg )
			if len(mol.entries)>0:
				self.mols[mol.mol_id] = mol

		# Generate map between ZMW and molecule IDs (self.zmw_mol_map)
		self.mols, self.zmw_mol_map = generate_molecule_ZMW_map( self.mols, self.chunk_id)

		# Exclude any molecules that are divided between split-up alignment files
		self.mols = remove_split_up_molecules( self.mols, self.split_mols )

		# Generate the IPD arrays per genomic position/strand by aggregating all 
		# IPD entries across molecules (self.ipdArrays)
		self.ipdArrays = self.create_agg_IPD_arrays()

		reader.close()

		# Return the processed IPD array dictionary
		return self.ipdArrays

	def create_agg_IPD_arrays( self ):
		"""
		Given a set of molecule object, each with a set of molecule identifiers
		and entries (containing IPD, base, pos, strand, call, subread), construct
		a dictionary of arrays that aggregate all the molecules together.
		"""
		tmp_ipdArrays  = {0:defaultdict(list), 1:defaultdict(list)}
		ipdArrays      = {0:{}, 1:{}}
		for mol in self.mols.values():
			for entry in mol.entries.values():
				if entry.call != "*":
					tmp_ipdArrays[entry.strand][entry.pos].append(entry.ipd)

		for strand in tmp_ipdArrays.keys():
			for pos in tmp_ipdArrays[strand].keys():
				if len(tmp_ipdArrays[strand][pos]) > 0:
					ipdArrays[strand][pos] = np.array(tmp_ipdArrays[strand][pos])
		return ipdArrays

class native_molecules_processor:
	def __init__( self,                   \
				  cmph5,                  \
				  chunk_id,               \
				  prefix,                 \
				  opts,                   \
				  fastq,                  \
				  ref,                    \
				  movie_name_ID_map,      \
				  control_ipds,           \
				  ref_size,               \
				  sites_pos,              \
				  sites_neg,              \
				  idx,                    \
				  split_mols):
		self.cmph5                   = cmph5
		self.chunk_id                = chunk_id
		self.prefix                  = prefix
		self.fastq                   = fastq
		self.ref                     = ref
		self.movie_name_ID_map       = movie_name_ID_map
		self.control_ipds            = control_ipds
		self.ref_size                = ref_size
		self.sites_pos_fn            = sites_pos
		self.sites_neg_fn            = sites_neg
		self.idx                     = idx
		self.opts                    = opts
		self.split_mols              = split_mols

		self.contig_id               = opts.contig_id
		self.align                   = opts.align
		self.nativeCovThresh         = opts.nativeCovThresh
		self.first_skip              = opts.firstBasesToSkip
		self.last_skip               = opts.lastBasesToSkip
		self.upstreamSkip            = opts.upstreamSkip
		self.downstreamSkip          = opts.downstreamSkip
		self.SMp                     = opts.SMp
		self.leftAnchor              = opts.leftAnchor
		self.rightAnchor             = opts.rightAnchor
		self.nat_lib                 = opts.nat_lib
		self.wgaCovThresh            = opts.wgaCovThresh
		self.out                     = opts.out

		if self.sites_pos_fn != None and self.sites_neg_fn != None:
			logging.debug("Process %s: reading motif sites..." % self.chunk_id)
			self.sites_pos = read_in_motif_sites( self.sites_pos_fn  )
			self.sites_neg = read_in_motif_sites( self.sites_neg_fn  )

	def __call__( self ):
		reader = CmpH5Reader(self.cmph5)

		self.chunk_output_fn = "%s_%s.tmp" % (self.out, self.chunk_id)
		self.var_chunk_fn    = "vars_%s.tmp" % self.chunk_id
		if self.align:
			self.var_f       = open(self.var_chunk_fn, "w")

		mol_alignments = defaultdict(list)
		for i,alignment in enumerate(reader[self.idx]):
			if alignment.accuracy >= self.opts.minAcc and alignment.MapQV >= self.opts.minMapQV and len(alignment.alignmentArray()) >= self.opts.minSubreadLen:
				mol_id = alignment.MoleculeID
				mol_alignments[mol_id].append(alignment)

		self.mols = {}
		incr      = int(max(4,math.floor(float(len(mol_alignments.keys())))/4))
		for i,(mol_id, alignments) in enumerate(mol_alignments.iteritems()):
			if i%incr==0:
				logging.info("...chunk %s - processing molecules: %s/%s (%.1f%%)" % (self.chunk_id, i, len(mol_alignments.keys()), 100*i/len(mol_alignments.keys())))
			mol = molecule( alignments, self.prefix, self.leftAnchor, self.rightAnchor, self.contig_id, self.sites_pos, self.sites_neg )
			if len(mol.entries)>0:
				self.mols[mol.mol_id] = mol

		# Exclude any molecules that are divided between split-up alignment files
		self.mols = remove_split_up_molecules( self.mols, self.split_mols )

		# [Optional]: align CCS reads to reference to find SNPs/errors
		if not self.align:
			# Need to empirically try to determine subread start/end positions in order to designate off-limits entries.
			for mol in self.mols.values():
				mol.var_pos = []
				self.empirical_get_start_end_pos( mol )
		elif len(self.mols.values()) > 0:
			CCS = CCS_aligner.mols_aligner( self.mols,                \
											self.fastq,               \
											self.ref,                 \
											self.movie_name_ID_map	, \
											self.align,               \
											self.chunk_id)
			# Output the called CCS read-level variants/errors to a chunk file
			for mol in self.mols.values():
				vars_str = ",".join(map(lambda x: str(x), mol.var_no_sc))
				self.var_f.write("%s %s\n" % (mol.mol_id, vars_str))

		if self.nat_lib=="short":
			# If the empirical start/end discovery showed a lack of positions with sufficient coverage, remove molecule
			del_me = [mol.mol_id for mol in self.mols.values() if mol.to_del]
			logging.debug("Process %s (chunk %s): deleting %s molecules due to too many positions with low coverage." % (self.chunk_id, \
																														  i, \
																														  len(del_me)))
			for mol_id in del_me:
				del self.mols[mol_id]

		if len(self.mols.values()) > 0:
			# Identify and remove positions to be excluded from further analysis
			tot_entries         = 0
			tot_entries_deleted = 0
			for mol in self.mols.values():
				entries_deleted, entries = self.remove_off_limits_positions( mol )
				tot_entries         += entries
				tot_entries_deleted += entries_deleted
			pct_deleted = float(tot_entries_deleted) / tot_entries * 100
			logging.debug("Process %s (chunk %s): deleted %s (%.1f%%) off-limits positions." % (self.chunk_id, \
																								 i, \
																								 tot_entries_deleted, \
																								 pct_deleted))

		# Generate the IPD arrays per genomic position/strand by aggregating all IPD entries across molecules (self.ipdArrays)
		logging.debug("Process %s: generating IPD arrays..." % self.chunk_id)
		for mol in self.mols.values():
			self.create_arrays( mol )

		if self.SMp:
			for mol in self.mols.values():
				mol.ipdArrays = self.condense_native_mol_motifs_into_one_pos( mol )

		# Now run the comparison test
		logging.debug("Process %s: running comparisons..." % self.chunk_id)
		for mol in self.mols.values():
			self.get_scores( mol )

		mols_w_results = len([mol for mol in self.mols.values() if len(mol.output)>0])
		logging.debug("Process %s (chunk %s): %s molecules generated comparison test output" % (self.chunk_id, i, mols_w_results))

		chunk_mols_with_output = []
		self.chunk_dirname = "chunk%s" % self.chunk_id
		if os.path.exists(self.chunk_dirname): shutil.rmtree(self.chunk_dirname)
		os.mkdir(self.chunk_dirname)
		for mol in self.mols.values():
			if len(mol.output) > 0:
				self.print_output( mol )
				chunk_mols_with_output.append( mol.mol_id )

		if mols_w_results > 0:
			self.concatenate_mol_results()

		shutil.rmtree(self.chunk_dirname)

		if self.align:
			self.var_f.close()
		reader.close()
		return self.chunk_output_fn

	def condense_native_mol_motifs_into_one_pos( self, mol ):
		"""
		If self.SMp==True (for long-insert phasing studies at low per-molecule coverage), 
		take all the IPDs from the specified motif and condense them into one position/strand. Specifically,
		move them all to the smallest genomic position of the motif on the molecule. The only IPD arrays present
		at this point should be at the motif sites.
		"""
		new_ipdArrays                 = {0:{}, 1:{}}
		mol.track_IPD_donor_positions = {0:defaultdict(list), 1:defaultdict(list)}
		for strand in mol.ipdArrays.keys():
			if len(mol.ipdArrays[strand].keys()) >= self.nativeCovThresh:
				smallest_pos = min(mol.ipdArrays[strand].keys())
				new_ipdArray = mol.ipdArrays[strand][smallest_pos]
				for pos in mol.ipdArrays[strand].keys():
					if pos != smallest_pos:
						x            = new_ipdArray
						y            = mol.ipdArrays[strand][pos]
						new_ipdArray = np.concatenate( [x, y] )

					mol.track_IPD_donor_positions[strand][smallest_pos].append(pos)
				
				new_ipdArrays[strand][smallest_pos] = new_ipdArray
		return new_ipdArrays

	def read_in_high_confidence_mod_sites( self, highConfCalls_fn ):
		"""
		Given a list of molecules/positions/strands, read in the sites where you want to collect
		IPD information.
		"""
		highConfs = []
		for line in open(highConfCalls_fn).xreadlines():
			# strand pos mol p nat_cov
			line    = line.strip()
			strand  =   int(line.split(" ")[0])
			pos     =   int(line.split(" ")[1])
			mol     =   int(line.split(" ")[2])
			WGA_bin =   int(line.split(" ")[13])
			highConfs.append( (mol, pos, strand, WGA_bin) )
		return highConfs 

	def apply_per_mol_coverage_filter( self ):
		"""
		We only want to keep molecules with sufficient per-molecule/strand coverage
		"""
		logging.debug("Process %s: %s molecules prior to coverage filtering." % (self.chunk_id, len(self.mols.values())))
		mols_dict_covFiltered = {}
		for mol_id, mol in self.mols.iteritems():
			subreads_strands = set([(entry.subread, entry.strand) for entry in mol.entries.values()])
			plus_subreads    = len([ subread[0] for subread in subreads_strands if subread[1]==0 ])
			minus_subreads   = len([ subread[0] for subread in subreads_strands if subread[1]==1 ])
			# This used to take the minimum of the two, but should take the max so that
			# more molecules survive this initial coverage filter. A second coverage
			# filtering will skip t-test analysis of distributions with < nativeCovThresh 
			# IPDs in the distribution.
			mol.cov          = max(plus_subreads, minus_subreads)

			if mol.cov >= self.nativeCovThresh or self.nat_lib=="long":
				mols_dict_covFiltered[mol_id] = mol

		if len(mols_dict_covFiltered.keys())==0:
			logging.debug("Process %s: no native molecules made it through the coverage filter!" % self.chunk_id)
		self.mols = mols_dict_covFiltered
		logging.debug("Process %s: %s molecules have coverage >= %.1f" % (self.chunk_id, len(self.mols), self.nativeCovThresh))

	def empirical_get_start_end_pos( self, mol ):
		"""
		When skipping the CCS alignment step, we need to try to determine reference-based
		coordinates of the start and end of each molecule template.
		"""
		cov_counter = Counter()
		for entry in mol.entries.values():
			cov_counter[entry.pos] += 1

		# valid_positions = [pos for (pos, cov) in cov_counter.iteritems() if cov >= (2*self.nativeCovThresh)]
		# Want to allow for some IPDs being thrown out due to +1:-1 filtering when pulling from cmp.h5 file.
		# That's why I'm multiplying by 1.5 instead of 2 (looking at both strands combined here)
		valid_positions = [pos for (pos, cov) in cov_counter.iteritems() if cov >= (1.5*self.nativeCovThresh)]
		if self.nat_lib=="long":
			valid_positions = [pos for (pos, cov) in cov_counter.iteritems()]
		
		if len(valid_positions) == 0:
			mol.to_del = True
		else:	
			mol.first_pos   = min(valid_positions)
			mol.last_pos    = max(valid_positions)

	def remove_off_limits_positions( self, mol ):
		"""Make a list of positions around the SNP to avoid during analysis.

		Also identify the earliest base position in the read. Will want to skip
		the first seven bases because of the adapter context not matching the
		reference context and introducing false kinetic signatures."""
		plus_off_limits  = []
		minus_off_limits = []
		for var in mol.var_pos:
			plus_off_limits  += range(var-self.upstreamSkip,   var+self.downstreamSkip+1)
			minus_off_limits += range(var-self.downstreamSkip, var+self.upstreamSkip+1)

		# Skip the first 'first_skip' positions
		plus_off_limits  += range(mol.first_pos, mol.first_pos + self.first_skip + 1)
		minus_off_limits += range(mol.last_pos,  mol.last_pos  - self.first_skip - 1, -1)

		# Skip the last 'last_skip' positions
		plus_off_limits  += range(mol.last_pos, mol.last_pos - self.last_skip + 1, -1)
		minus_off_limits += range(mol.last_pos, mol.last_pos + self.last_skip - 1)

		plus_off_limits  = set(plus_off_limits)
		minus_off_limits = set(minus_off_limits)

		for entry in mol.entries.values():
			entry.off_limits = False
			if entry.strand == 0 and entry.pos in plus_off_limits:
				entry.off_limits = True
			elif entry.strand == 1 and entry.pos in minus_off_limits:
				entry.off_limits = True

		keys_to_del = []
		for key, entry in mol.entries.iteritems():
			if entry.off_limits:
				keys_to_del.append(key)

		before = len(mol.entries)
		for key in keys_to_del:
			del mol.entries[key]
		after  = len(mol.entries)
		return (before - after), before

	def create_arrays( self, mol ):
		"""
		Given a set of molecule object, each with a set of molecule identifiers
		and entries (containing IPD, base, pos, strand, call, subread), construct
		a dictionary of arrays that aggregate all the molecules together.
		"""
		ipdArrays                = {0:defaultdict(list), 1:defaultdict(list)}
		mol.ipdArray_subread_map = {}
		
		for entry in mol.entries.values():
			if entry.call != "*":
				position_in_array = len(ipdArrays[entry.strand][entry.pos])
				ipdArrays[entry.strand][entry.pos].append(entry.ipd)
				mol.ipdArray_subread_map[ (entry.strand, entry.pos, position_in_array) ] = entry.subread

		mol.ipdArrays = {0:{}, 1:{}}
		for strand in ipdArrays.keys():
			for pos in ipdArrays[strand].keys():
				if len(ipdArrays[strand][pos]) >= self.nativeCovThresh or self.nat_lib=="long":
					mol.ipdArrays[strand][pos] = np.array(ipdArrays[strand][pos])

	def get_scores( self, mol ):
		""""""
		mol.output   = []
		for strand in mol.ipdArrays.keys():
			if len(mol.ipdArrays[strand].keys()) > 0:
				for pos, array in mol.ipdArrays[strand].iteritems():
					if strand in set(self.control_ipds.keys()) and pos in set(self.control_ipds[strand].keys()):
						if len(array) < self.nativeCovThresh:
							# Native per-molecule coverage at this position/strand is too low -- skip.
							continue

						if self.SMp:
							# Perform the same combination of IPDs as we did on the native molecule
							positions_to_combine = mol.track_IPD_donor_positions[strand][pos]
							control_ipds         = self.control_ipds[strand][pos]
							for combine_pos in positions_to_combine:
								try:
									to_add       = self.control_ipds[strand][combine_pos]
									control_ipds = np.concatenate( [control_ipds, to_add] )
								except KeyError:
									# IPDs were pooled from this position in the native molecule, but no
									# control IPDs are present here to add to the control pool
									continue

						else:
							control_ipds = self.control_ipds[strand][pos]

						control_ln_vals  = map(lambda x: math.log(x + 0.0001), control_ipds)

						if len(control_ipds) < self.wgaCovThresh:
							# WGA coverage at that location is too low -- skip.
							continue
					
						# Convert IPD values to ln(IPD) to get into normal distribution
						control_IPD_mean = np.array(control_ln_vals).mean()

						array = np.array(map(lambda x: math.log(x+0.0001), array))
						
						line = "%s\t%s\t%.3f\t%s\t%.3f\t%.3f\t%s\t%s\t%s" % (strand,            \
															   pos,               \
															   array.mean() - control_IPD_mean, \
															   mol.mol_id,        \
															   array.mean(),      \
															   control_IPD_mean,  \
															   len(array),        \
															   len(control_ln_vals), \
															   np.abs(mol.last_pos - mol.first_pos))
						mol.output.append(line)

					else:
						logging.debug("No entry for position %s in local control dict)" % pos)

	def print_output( self, mol ):
		"""
		Print the output lines to a molecule-specific tmp file. These will be catted together
		at the end.
		"""
		out_fn  = "%s/mol_%s.tmp.txt" % (self.chunk_dirname, mol.mol_id)
		outfile = open(out_fn, "w")
		for line in mol.output:
			outfile.write("%s\n" % line)
		outfile.close()

	def concatenate_mol_results( self ):
		"""
		"""
		fns       = glob.glob("%s/mol_*.tmp.txt" % self.chunk_dirname)
		# fns       = map(lambda x: "mol_%s.tmp.txt" % x, mols_with_output)
		cat_CMD   = "cat %s >> %s" % (" ".join(fns), self.chunk_output_fn)
		p         = subprocess.Popen(cat_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdOutErr = p.communicate()
		sts       = p.returncode
		if sts != 0:
			raise Exception("Failed alignment command: %s" % cat_CMD)

		for fn in fns:
			try:
				os.remove(fn)
			except OSError:
				pass
