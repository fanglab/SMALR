import sys,os
import shutil
import numpy as np
import logging
import subprocess
import sam_parser

class mols_aligner:
	"""
	This class will coordinate the execution of the CCS read alignments.
	"""
	def __init__( self, mols, orig_fastq, ref, nat_movie_name_ID_map, align, chunk_id ):
		self.orig_fastq = orig_fastq
		self.ref        = ref
		self.chunk_id   = chunk_id

		sam_md, tmp_dir = self.align_CCS_mols_to_ref( mols.values(), nat_movie_name_ID_map, align )
		self.parse_alignments( sam_md, mols, align )
		# self.native_stats = self.get_native_read_stats( mols )
		shutil.rmtree(tmp_dir)

	def align_CCS_mols_to_ref( self, mols, nat_movie_name_ID_map, align ):
		"""Pull out the CCS reads for every molecule of interest and generate on fastq file.
		Then align this file to the reference using bwa bwasw. Finally, transform the output
		so that it returns a SAM file with the necessary MD field (for mismatch calling)."""
		self.readname_molID_map = {}

		def get_CCS_reads_for_aligning( mols, nat_movie_name_ID_map, align ):
			"""Go through the CCS fastq and index where each molecule/ZMW resides. Should speed
			up access when pulling out the reads."""
			tmp_dir = "tmp_align_%s" % self.chunk_id
			if os.path.exists(tmp_dir): 
				shutil.rmtree(tmp_dir)
			os.mkdir(tmp_dir)

			chunk_fastq = "chunk%s_mols_CCS.fastq" % self.chunk_id
			new_fastq   = os.path.join(tmp_dir, chunk_fastq)
			new_fastq_f = open(new_fastq,"w")

			id_fn     = os.path.join(tmp_dir, "all_ccs.fastq.index")
			CMD       = "grep -n ^@m %s > %s" % (self.orig_fastq, id_fn)
			p         = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdOutErr = p.communicate()
			sts       = p.returncode
			if os.path.exists(id_fn) and os.path.getsize(id_fn) > 0:
				zmw_line_map     = {}
				zmw_readname_map = {}
				for line in open(id_fn).xreadlines():
					line       = line.strip()
					line_n     = int(line.split(":")[0])
					readname   =     line.split(":")[1].strip("@")
					movie_name = readname.split("/")[0]

					try:
						movie_id = nat_movie_name_ID_map[movie_name]
					except KeyError:
						# The queried read is not in this chunk of molecules
						continue
					if readname.split("/")[-1] == "ccs":
						zmw_id                   = readname.split("/")[-2] + "_" + str(movie_id)
					else:
						zmw_id                   = readname.split("/")[-1] + "_" + str(movie_id)
					zmw_line_map[zmw_id]     = line_n
					zmw_readname_map[zmw_id] = readname

				n_mols       = len(mols)
				not_in_fastq = 0
				for n,mol in enumerate(mols):
					try:
						target    = zmw_line_map[mol.zmw_id]
						i = 0
						with open(self.orig_fastq) as f:
							
							while i < target-1:
								f.next()
								i += 1
							j = 0
							while j < 4:
								# print f.next()
								new_fastq_f.write(f.next())
								j += 1
							mol.in_fastq = True
							readname     = zmw_readname_map[mol.zmw_id]
							self.readname_molID_map[readname] = mol.mol_id
					except KeyError:
						not_in_fastq += 1
						progress      = float(n)/n_mols * 100
						logging.debug("Chunk %s: ZMW %s (mol %s) not in fastq file. Progress -- %.1f%%" % (self.chunk_id, mol.zmw_id, mol.mol_id, progress))
						pass
				pct_not_found = float(not_in_fastq) / len(mols) * 100
				logging.debug("Chunk %s: %s reads not found (%.1f%%) in original CCS fastq file." % (self.chunk_id, not_in_fastq, pct_not_found))
			else:
				logging.error("CCS fastq indexing failed: %s" % CMD)
			new_fastq_f.close()
			return new_fastq, tmp_dir

		def call_bwasw_samtools( fq ):
			
			def run_command( CMD ):
				p         = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				stdOutErr = p.communicate()
				sts       = p.returncode
				if sts != 0:
					for entry in stdOutErr:
						print entry
					raise Exception("Failed alignment command: %s" % CMD)


			sam        = fq+".aligned.sam"
			bam        = fq+".aligned.bam"
			sam_md     = fq+".aligned.md.sam"

			bwa_CMD     = "bwa bwasw %s %s > %s" % (self.ref, fq, sam)
			run_command( bwa_CMD )

			# Check if there are any alignments in the file
			if os.path.getsize(sam) == 0:
				logging.error("No alignments for %s!" % fq)
				os.remove(sam)
			else:
				samtools_view_CMD  = "samtools view -b -S %s > %s" % (sam, bam)
				run_command(samtools_view_CMD)

				# We are going to use samtools calmd to calculate the MD field, which along
				# with the CIGAR string in the sam file will let us identify SNPs/Indels.
				samtools_calmd_CMD = "samtools calmd %s %s > %s" % (bam, self.ref, sam_md)
				run_command(samtools_calmd_CMD)
			return sam_md

		new_fastq, tmp_dir = get_CCS_reads_for_aligning( mols, nat_movie_name_ID_map, align )
		sam_md = None
		if os.path.getsize(new_fastq) > 0:
			sam_md    		   = call_bwasw_samtools( new_fastq )
		return sam_md, tmp_dir

	def parse_alignments( self, sam_md, mols, align ):
		"""Parse each line in the SAM file and add certain values from the alignment as
		attributes of the molecule object."""
		i = 0
		j = 0
		to_del = []
		if sam_md != None:
			for line in open(sam_md).xreadlines():
				if line[0]=="@": # header
					continue
				readname    = line.split("\t")[0]
				mol_id      = self.readname_molID_map[readname]
				if line.split("\t")[2] == "*": # No alignment found
					to_del.append(mol_id)
					continue

				# Do a quick check to make sure were are looking that the alignment
				# that corresponds to the ref positions listed in the IPD csv file.
				if readname.find("/") < 0:
					# This is a weird alignment where the readname is corrupt.
					j += 1
					to_del.append(mol_id)
					continue
				align_start = int(line.split("\t")[3])
				one_pos     = mols[mol_id].entries.values()[0].pos
				if np.abs(one_pos - align_start) > 2000:
					# We probably aren't looking at the proper alignment
					# Even if we are, there's a ton of soft-clipping at the
					# beginning of this read -- suspicious.
					logging.debug("Chunk %s: suspicious positions! mol:%s align_start:%s sampled_pos:%s" % (self.chunk_id, mol_id, align_start, one_pos))
					i += 1
					to_del.append(mol_id)
					continue

				readname, var_pos, softclip, CIGAR, MD, align_start, align_end = sam_parser.read_sam_line(line)
				if not align:
					var_pos = []
				mol_id          = self.readname_molID_map[readname]
				mol             = mols[mol_id]
				mol.mapped      = True
				mol.var_pos     = var_pos
				mol.softclip    = softclip
				mol.var_no_sc   = list(set(mol.var_pos) - set(mol.softclip))
				mol.CIGAR       = CIGAR
				mol.MD          = MD
				mol.align_start = align_start
				mol.align_end   = align_end
				mol.first_pos   = align_start
				mol.last_pos    = align_end
		logging.debug("Chunk %s: skipped %s alignments -- start pos far away from sampled position." % (self.chunk_id, i))
		logging.debug("Chunk %s: skipped %s alignments -- readname was just the movie name (?)." % (self.chunk_id, j))

		to_del = to_del + [mol.mol_id for mol in mols.values() if (not mol.mapped or not mol.in_fastq)]
		for mol_id in to_del:
			try:
				del mols[mol_id]
			except KeyError:
				pass

	def get_native_read_stats( self, mols ):
		"""
		"""
		all_mols           = []
		in_fastq_no_mapped = []
		no_mapped          = []
		no_fastq           = []
		to_analyze         = []
		for mol in mols.values():
			all_mols.append(mol.mol_id)
			if not mol.mapped:
				no_mapped.append(mol.mol_id)
			if not mol.in_fastq:
				no_fastq.append(mol.mol_id)
			if mol.mapped and mol.in_fastq:
				to_analyze.append(mol.mol_id)
		native_stats = {"all_mols"       : all_mols,      \
						"no_mapped"          : no_mapped, \
						"no_fastq"           : no_fastq,  \
						"to_analyze"         : to_analyze}
		return native_stats
