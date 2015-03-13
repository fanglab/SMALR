# Use this as a sequence to map an encoded operation to the appropriate
# character
DECODE = 'MIDNSHPC'

# This dictionary maps operations to their integer encodings
_ENCODE = dict( (c,i) for (i, c) in enumerate(DECODE) )

def parse_cigar( cigar ):
	"""
	Parse CIGAR string and return a list of pysam-encoded tuples.

	MIDNSHP => 0123456

	>>> parse("3S17M8D4M9I3H")
	[(4, 3), (0, 17), (2, 8), (0, 4), (1, 9), (5, 3)]
	"""
	result = []
	n = ''
	for c in cigar:
		if c.isdigit():
			n += c
		elif c in _ENCODE:
			if n == '':
				raise ValueError("end of CIGAR string reached, but an operator was expected")
			result.append( (_ENCODE[c], int(n)) )
			n = ''
	return result

def parse_md( md ):
	"""
	Parse the MD field to identify the mismatches in the SAM alignment.
	It looks like MD fields can also spot deletions, but not insertions.
	Luckily, we can get both insertions and deletions from the CIGAR.
	
	39^A98^C40C36^C22^C10^C56^G79^C44

	"""
	result = []
	md     = md.split(":")[-1]
	n      = ""
	DEL    = False
	MM     = False
	for i,c in enumerate(md):
		if c.isdigit():
			if DEL:
				result.append( (_ENCODE["D"], len(del_bases)) )
				n = ""
			if MM:
				result.append( (_ENCODE["C"], len(mm_bases)) )
				n = ""
			# Generating matching interval number
			DEL = False
			MM  = False
			n  += c
			mm_bases = []
			if i == (len(md)-1):
				# If you've reached the end of the MD string
				result.append( (_ENCODE["M"], int(n)) )
		elif c == "^":
			result.append( (_ENCODE["M"], int(n)) )
			DEL       = True
			del_bases = []
		elif c in ["A","C","G","T"]:
			if DEL:
				# Deleted sequence
				del_bases.append(c)
			else:
				# Mismatch
				if not MM:
					result.append( (_ENCODE["M"], int(n)) )
				MM = True
				mm_bases.append(c)
	# print "MD result", result
	return result

def read_sam_line( line ):
	"""
	We can extract all the positions of errors in the CCS alignment by
	querying both the CIGAR string (for indels) and the MD field (for 
	mismatches and deletions).

	>ref000001|NC_012920_HUMAN_MTDNA_NON
	GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGG
	GTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTC

	(1) query:
	@TEST_SEQ
	GTATGCACGTGATAGCATTCGAGACGCTGAAGCCGGAGCTACCCTATGTCGCAGTATCTGTCTTTGATTC
	         ^SNP      ^DEL      ^SNP      ^INS

	CIGAR = 19M1D20M1I30M
	MD    = MD:Z:9C9^G10G39
	"""
	# constants (including 7th 'C' coding for mismatches found in the MD field)
	M = 0; I = 1; D = 2; N = 3; S = 4; H = 5; P = 6; C = 7;

	var_pos  = []


	readname =     line.split("\t")[0]
	start    = int(line.split("\t")[3])
	QV       = int(line.split("\t")[4])
	CIGAR    =     line.split("\t")[5]
	MD       =     line.split("\t")[-1].strip()
	
	cig_result  = parse_cigar(CIGAR)
	pos         = start
	last_pos    = start
	softclip    = []
	for i,x in enumerate(cig_result):
		char = DECODE[x[0]]
		# Advance in ref coords by the length of this stretch
		pos += x[1]
		if   char == "M":
			pass
		elif char == "I":
			var_pos += range(last_pos, pos)
		elif char == "D":
			var_pos += range(last_pos, pos)
		elif char == "S":
			if i == 0:
				# This is the first element encountered.
				# The soft-clipped beginning sequence of the read is before the alignment
				# technically starts and should be marked as a var_pos
				pos = start
				softclip += range(start-1, start-x[1]-1, -1)
				var_pos  += range(start-1, start-x[1]-1, -1)
			else:
				softclip += range(last_pos, pos)
				var_pos  += range(last_pos, pos)
		else:
			raise Exception("Strange character in CIGAR %s: %s" % (CIGAR, char))
		last_pos = pos

	md_result = parse_md( MD )
	pos       = start
	last_pos  = start
	for x in md_result:
		char = DECODE[x[0]]
		# Advance in ref coords by the length of this stretch
		pos += x[1]
		if   char == "M":
			pass
		elif char == "C":
			# var_pos.append(pos)
			var_pos += range(last_pos, pos)
		elif char == "D":
			# var_pos.append(pos)
			var_pos += range(last_pos, pos)
		else:
			raise Exception("Strange character in MD %s: %s" % (MD, char))
		last_pos = pos

	return readname, list(set(var_pos)), list(set(softclip)), CIGAR, MD, start, last_pos