#! /usr/bin/python
###
###
###
### Script written by Ivan Sovic, August 2015.
### All rights reserved.
###
###
###
# Some comments on DALIGNER.
#
# Important about DB generation: DALIGNER splits the input reference fasta into chunks of non 'N' bases, as well as into chromosomes.
# E.g. Let's say we have hg19_chr3 as a reference. That is only one chromosome (or, only one sequence). However, fasta2DAM will generate 21 reference sequence.
# This corresponds to the number of breakpoints present in the chr3.
# As another example, if you convert S. Cerevisiae S288c, you will get exactly 16 sequences as expected.



import re;

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../src');

import subprocess;
import multiprocessing;

try:
	import basicdefines;
	USE_BASICDEFINES_ = True;
	ALIGNERS_PATH_ROOT_ABS_ = basicdefines.ALIGNERS_PATH_ROOT_ABS;
except:
	USE_BASICDEFINES_ = False;
	ALIGNERS_PATH_ROOT_ABS_ = SCRIPT_PATH;

ALIGNER_URL = 'https://github.com/thegenemyers/DALIGNER.git';
ALIGNER_DB_URL = 'https://github.com/thegenemyers/DAZZ_DB.git';
ALIGNER_PATH = ALIGNERS_PATH_ROOT_ABS_ + '/DALIGNER/DALIGNER/';
ALIGNER_DB_PATH = ALIGNERS_PATH_ROOT_ABS_ + '/DALIGNER/DAZZ_DB/';
BIN = 'daligner';
MAPPER_NAME = 'DALIGNER';

RUNNING_PATH = os.path.dirname(sys.argv[0]);

def peek(fp, num_chars):
	data = fp.read(num_chars);
	if len(data) == 0:
		return '';
	fp.seek(num_chars * -1, 1);
	return data;

# Returns a single read from the given FASTA/FASTQ file.
# Parameter header contains only the header of the read.
# Parameter lines contains all lines of the read, which include:
# - header
# - seq
# - '+' if FASTQ
# - quals if FASTQ
# Parameter lines is an array of strings, each for one component.
# Please note that multiline FASTA/FASTQ entries (e.g. sequence line)
# will be truncated into one single line.
def get_single_read(fp):
	lines = [];
	
	line = fp.readline();
	header = line.rstrip();
	header_leading_char = '';
	if (len(header) > 0):
		sequence_separator = header[0];
		header_leading_char = header[0];
		header = header[1:];			# Strip the '>' or '@' sign from the beginning.
	else:
		return ['', []];
	
	next_char = peek(fp, 1);
	
	line_string = '';
	lines.append(header_leading_char + header);
	
	num_lines = 1;
	#while len(next_char) > 0 and next_char != sequence_separator or (next_char == '@' and num_lines < 4):
	while (len(next_char) > 0 and (next_char != sequence_separator or (next_char == '@' and num_lines < 4))):
		line = fp.readline();
		if (line.rstrip() == '+' or line.rstrip() == ('+' + header)):
		#if (line.rstrip()[0] == '+'):
			lines.append(line_string);
			lines.append(line.rstrip());
			line_string = '';
		else:
			line_string += line.rstrip();
		next_char = peek(fp, 1);
		num_lines += 1;
		
	lines.append(line_string);
	
	return [header, lines];

def get_fastq_headers_and_lengths(fastq_path):
	headers = [];
	lengths = [];
	### DALIGNER splits each read/sequence into 'subreads' at every N position. That's why bread index does not necessarily correspond to the actual sequence ID.
	daligner_seq_id = [];	### Contains an array of tuples (seq_id, seq_header, start_offset).

	fp_in = None;

	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % fastq_path;
		exit(1);
	
	seq_id = 0;

	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		headers.append(header);
		lengths.append(len(read[1]));

		seq = read[1];

		i = 0;
		while (i < len(seq)):
			if (seq[i] != 'N' and (i == 0 or (i > 0 and seq[i-1] == 'N'))):
				daligner_seq_id.append( (seq_id, header, i) );
			i += 1;
		# i = 0;
		# while (i < len(seq)):
		# 	if (i == 0 or (i > 0 and seq[i] == 'N' and seq[i-1] != 'N')):
		# 		daligner_seq_id.append( (seq_id, header, i) );
		# 	i += 1;
		
		seq_id += 1;

	fp_in.close();

	# print 'len(daligner_seq_id) = %d' % (len(daligner_seq_id));
	# exit(1);

	return [headers, lengths, daligner_seq_id];



def measure_command_wrapper(out_filename):
	if (USE_BASICDEFINES_ == True):
		return basicdefines.measure_command(out_filename);
	else:
		return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % out_filename;

def parse_memtime(memtime_path):
	cmdline = '';
	realtime = 0;
	cputime = 0;
	usertime = 0;
	systemtime = 0;
	maxrss = 0;
	rsscache = 0;
	time_unit = '';
	mem_unit = '';

	try:
		fp = open(memtime_path, 'r');
		lines = [line.strip() for line in fp.readlines() if (len(line.strip()) > 0)];
		fp.close();
	except Exception, e:
		sys.stderr.write('Could not find memory and time statistics in file "%s".\n' % (memtime_path));
		return [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit];

	for line in lines:
		if (line.startswith('Command line:')):
			cmdline = line.split(':')[1].strip();
		elif (line.startswith('Real time:')):
			split_line = line.split(':')[1].strip().split(' ');
			realtime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('CPU time:')):
			split_line = line.split(':')[1].strip().split(' ');
			cputime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('User time:')):
			split_line = line.split(':')[1].strip().split(' ');
			usertime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('System time:')):
			split_line = line.split(':')[1].strip().split(' ');
			systemtime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('Maximum RSS:')):
			split_line = line.split(':')[1].strip().split(' ');
			maxrss = float(split_line[0].strip());
			mem_unit = split_line[1].strip();
		# elif (line.startswith('')):
		# 	split_line = line.split(':')[1].strip().split(' ');
		# 	rsscache = float(split_line[0].strip());
		# 	mem_unit = split_line[1].strip();

	return [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit];

def parse_memtime_files_and_accumulate(memtime_files, final_memtime_file):
	final_command_line = '';
	final_real_time = 0.0;
	final_cpu_time = 0.0;
	final_user_time = 0.0;
	final_system_time = 0.0;
	final_time_unit = '';
	final_max_rss = 0;
	final_mem_unit = '';

	i = 0;
	for memtime_file in memtime_files:
		i += 1;
		sys.stderr.write('Parsing memtime file "%s"...\n' % (memtime_file));

		[cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit] = parse_memtime(memtime_file);
		if (i == 1):
			final_command_line = cmdline;
			final_real_time = realtime;
			final_cpu_time = cputime;
			final_user_time = usertime;
			final_system_time = systemtime;
			final_max_rss += maxrss;
			final_time_unit = time_unit;
			final_mem_unit = mem_unit;
		else:
			if (time_unit == final_time_unit and mem_unit == final_mem_unit):
				final_command_line += '; ' + cmdline;
				final_real_time += realtime;
				final_cpu_time += cputime;
				final_user_time += usertime;
				final_system_time += systemtime;
				final_max_rss += maxrss;
			else:
				sys.stderr.write('Memory or time units not the same in all files! Instead of handling this, we decided to be lazy and just give up.\n');
				break;

	try:
		fp = open(final_memtime_file, 'w');
	except Exception, e:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % (final_memtime_file));
		return;

	if (final_cpu_time <= 0.0):
		final_cpu_time = final_user_time + final_system_time;

	fp.write('Command line: %s\n' % (final_command_line));
	fp.write('Real time: %f %s\n' % (final_real_time, final_time_unit));
	fp.write('CPU time: %f %s\n' % (final_cpu_time, final_time_unit));
	fp.write('User time: %f %s\n' % (final_user_time, final_time_unit));
	fp.write('System time: %f %s\n' % (final_system_time, final_time_unit));
	fp.write('Maximum RSS: %f %s\n' % (final_max_rss, final_mem_unit));

	fp.close();





# Command line: /home/isovic/work/eclipse-workspace/git/aligneval/src/../tools/cgmemtime/cgmemtime -o /home/isovic/work/eclipse-workspace/git/aligneval/src/../evaluation/reads-simulated/OxfordNanopore-pbsim-observed_graphmap-1d-1k/saccharomyces_cerevisiae/GraphMap-v1.memtime /home/isovic/work/eclipse-workspace/git/aligneval/wrappers/../aligners/graphmap/bin/Linux-x64//graphmap -v 5 -t 4 -B 0 -b 3 -r /home/isovic/work/eclipse-workspace/git/aligneval/src/../reference-genomes/saccharomyces_cerevisiae.fa -d /home/isovic/work/eclipse-workspace/git/aligneval/src/../reads-simulated/OxfordNanopore-pbsim-observed_graphmap-1d-1k/saccharomyces_cerevisiae/reads.fa -o /home/isovic/work/eclipse-workspace/git/aligneval/src/../evaluation/reads-simulated/OxfordNanopore-pbsim-observed_graphmap-1d-1k/saccharomyces_cerevisiae/GraphMap-v1.sam
# Real time:   18.164 s
# CPU time:   61.911 s
# User time:   54.086 s
# System time:    7.825 s
# Maximum RSS: 1327 MB

def read_fastq(fastq_path):
	headers = [];
	seqs = [];
	quals = [];
	
	fp_in = None;

	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % fastq_path;
		return;
	
	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		
		seq = read[1];
		qual = '';
		if (len(read) == 4):
			qual = read[3];
		headers.append(header);
		seqs.append(seq);
		quals.append(qual);
		
	fp_in.close();
	
	return [headers, seqs, quals];

def convert_to_fasta(fastq_path, out_fasta_path):
	headers = [];
	seqs = [];
	quals = [];
	fp_in = None;
	fp_out = None;
	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % fastq_path;
		return;
	try:
		fp_out = open(out_fasta_path, 'w');
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % out_fasta_path;
		fp_in.close();
		return;
	while True:
		[header, read] = get_single_read(fp_in);
		if (len(header) == 0):
			break;
		seq = read[1];
		fp_out.write('>' + header + '\n');
		fp_out.write(seq + '\n');
	fp_in.close();
	fp_out.close();



def wrap_fasta_file(fasta_file, daligner_fasta_file):
	try:
		fp_in = open(fasta_file, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % fasta_file);
		exit(0);

	try:
		fp_out = open(daligner_fasta_file, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % daligner_fasta_file);
		exit(0);

	current_read = 0;

	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		current_read += 1;

		if (len(read[1]) <= 20):	### DALIGNER has a lower length limit of 10bp.
			continue;

		read[1] = re.sub("(.{500})", "\\1\n", read[1], 0, re.DOTALL);	### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
		if (len(read) == 4):
			read[3] = re.sub("(.{500})", "\\1\n", read[3], 0, re.DOTALL);	### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();
	fp_out.close();

def convert_reads_to_pacbio_format(reads_file, daligner_reads_file):
	try:
		fp_in = open(reads_file, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % reads_file);
		exit(0);

	try:
		fp_out = open(daligner_reads_file, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % daligner_reads_file);
		exit(0);

	current_read = 0;

	header_conversion_hash = {};

	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		current_read += 1;

		if (len(read[1]) <= 10):	### DALIGNER has a lower length limit of 10bp.
			sys.stderr.write('Found a read shorter than 10bp. Removing from the output.\n');
			continue;

		### Check if the read is already formatted like PacBio.
		if (header.count('/') == 2 and 'RQ' in header):
			fp_out.write('\n'.join(read) + '\n');
			continue;

		trimmed_header = header.replace('_', ' ').split()[0];
		# pacbio_header = '%s/%d/0_%d RQ=0.850' % (trimmed_header, current_read, len(read[1]));
		pacbio_header = 'S1/%d/0_%d RQ=0.850' % (current_read, len(read[1]));
		header_conversion_hash[pacbio_header] = header;
		read[0] = '%s%s' % (read[0][0], pacbio_header); ### Keep the first char of the header line.
		read[1] = re.sub("(.{500})", "\\1\n", read[1], 0, re.DOTALL);	### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
		if (len(read) == 4):
			read[3] = re.sub("(.{500})", "\\1\n", read[3], 0, re.DOTALL);	### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();
	fp_out.close();

	return header_conversion_hash;

def execute_command(command):
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_get_stdout(command):
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	[out, err] = p.communicate()
	sys.stderr.write('\n');

	return [out, err];

class Overlap:
	def __init__(self, original_line=None, match_obj=None):
		self.bread = 0;		# For alignments this is the ID of the reference hit. It's 1-based.
		self.aread = 0;		# ID of the read that is aligned/overlapped. It's 1-based.
		self.orient = 0;	# 'n' or 'c'.
		self.bstart = 0;	# B-read/reference start.
		self.bend = 0;		# B-read/reference end.
		self.astart = 0;	# A-read start.
		self.aend = 0;		# A-read end.
		self.diffs = 0;		# Number of diffs.
		self.tracepts = 0;	# Number of trace points.
		self.aln_ref = '';
		self.aln_query = '';
		self.aln_matching = '';
		self.original_line = '';

		if (match_obj != None):
			self.assign(match_obj, original_line);

	def assign(self, match_obj, original_line=None):
		[self.bread, self.aread, self.orient, self.bstart, self.bend, self.astart, self.aend, self.diffs, self.tracepts] = match_obj.groups();
		self.bread = int(''.join(self.bread.split(',')));
		self.aread = int(''.join(self.aread.split(',')));
		self.bstart = int(''.join(self.bstart.split(',')));
		self.bend = int(''.join(self.bend.split(',')));
		self.astart = int(''.join(self.astart.split(',')));
		self.aend = int(''.join(self.aend.split(',')));
		self.diffs = int(''.join(self.diffs.split(',')));
		self.tracepts = int(''.join(self.tracepts.split(',')));
		self.aln_ref = '';
		self.aln_query = '';
		self.aln_matching = '';
		self.original_line = '';
		if (original_line != None):
			self.original_line = original_line;

	def verbose_as_string(self):
		ret = '';
		ret += 'bread = %d\n' % (self.bread);						### Reference hit ID
		ret += 'aread = %d\n' % (self.aread);						### Read ID
		ret += 'orient = %s\n' % (self.orient);						### 'n' or 'c'
		ret += 'bstart = %d\n' % (self.bstart);
		ret += 'bend = %d\n' % (self.bend);
		ret += 'astart = %d\n' % (self.astart);
		ret += 'aend = %d\n' % (self.aend);
		ret += 'diffs = %d\n' % (self.diffs);
		ret += 'tracepts = %d\n' % (self.tracepts);
		ret += 'aln_ref =      %s\n' % (self.aln_ref);
		ret += 'aln_matching = %s\n' % (self.aln_matching);
		ret += 'aln_query =    %s\n' % (self.aln_query);
		return ret;

	def add_ref_alignment(self, ref_aln):
		self.aln_ref += ref_aln + '\n';

	def add_matching_alignment(self, matching_aln):
		self.aln_matching += matching_aln + '\n';

	def add_query_alignment(self, query_aln):
		self.aln_query += query_aln + '\n';

	def calc_cigar_string(self, read_length, cigar_format='basic'):
		if (cigar_format == 'basic'):
			TB_MATCH = 'M'; # 0;
			TB_MISMATCH = 'M'; # 1;
			TB_INSERTION = 'I'; # 2;
			TB_DELETION = 'D'; # 3;
		else:
			TB_MATCH = '='; # 0;
			TB_MISMATCH = 'X'; # 1;
			TB_INSERTION = 'I'; # 2;
			TB_DELETION = 'D'; # 3;

		traceback = [];

		### Convert the visual alignment to a traceback array.
		i = 0;
		while (i < len(self.aln_matching)):
			if (self.aln_matching[i] == '|'):
				traceback.append(TB_MATCH);
			elif (self.aln_matching[i] == '*'):
				if (self.aln_ref[i] == '-' and self.aln_query[i] != '-'):
					traceback.append(TB_INSERTION);
				elif (self.aln_ref[i] != '-' and self.aln_query[i] == '-'):
					traceback.append(TB_DELETION);
				elif (self.aln_ref[i] != self.aln_query[i] and self.aln_ref[i] != '-' and self.aln_query[i] != '-'):
					traceback.append(TB_MISMATCH);
			i += 1;

		### Summarize the traceback array into a CIGAR string.
		cigar = [];
		last = ''
		results = []
		for op in traceback:
			if op == last:
			    cigar[-1] = (op, cigar[-1][1] + 1);
			else:
			    cigar.append((op, 1));
			    last = op;

		cigar_string = '%dS' % (self.astart) if (self.astart > 0) else '';
		cigar_string += ''.join(['%d%s' % (op[1], op[0]) for op in cigar]);
		cigar_string += '%dS' % (read_length - self.aend) if (read_length > (self.aend + 1)) else '';

		return cigar_string;



	def convert_to_sam(self, ref_headers, ref_daligner_seq_id, read_headers, read_seqs, read_quals, header_conversion_hash):
		# num_clip_front = int(qstart) - 1;
		# num_clip_back = int(qlen) - (int(qend));
		# sam_cigar = convert_btop_to_cigar(btop, num_clip_front, num_clip_back, sstrand);

		try:
			qname = header_conversion_hash[read_headers[self.aread-1]];
			# .split()[0];
		except:
			sys.stderr.write('ERROR: Read "%s" cannot be found in the reads file! Faulty alignment in DALIGNER output?\n');
			exit(1);

		# qname = header_conversion_hash[qname];

		daligner_seq_id = ref_daligner_seq_id[self.bread-1];	### Contains an array of tuples (seq_id, seq_header, start_offset).

		flag = 0 if (self.orient == 'n') else 16;
		rname = daligner_seq_id[1];
		pos = daligner_seq_id[2] + self.bstart + 1;
		mapq = 255;
		sam_seq = read_seqs[self.aread-1] if (self.orient == 'n') else revcomp_seq(read_seqs[self.aread-1]);
		sam_qual = read_quals[self.aread-1] if (self.orient == 'n') else read_quals[self.aread-1][::-1];
		sam_cigar = self.calc_cigar_string(len(sam_seq));
		if (len(sam_qual) == 0):
			sam_qual = '*';
		sam_NM = self.diffs;

		sam_line = '';
		sam_line += '%s\t' % (qname);						# 1. qname
		sam_line += '%d\t' % (flag);						# 2. flag
		sam_line += '%s\t' % (rname);						# 3. rname
		sam_line += '%d\t' % (pos);				# 4. pos
		sam_line += '%d\t' % (mapq);								# 5. mapq
		sam_line += '%s\t' % (sam_cigar);					# 6. CIGAR
		sam_line += '*\t';									# 7. rnext
		sam_line += '0\t';									# 8. pnext
		sam_line += '0\t';									# 9. tlen
		sam_line += '%s\t' % (sam_seq);						# 10. seq
		sam_line += '%s\t' % (sam_qual);					# 11. qual
		sam_line += 'AS:i:%d\t' % (self.aend - self.astart + 1 - sam_NM);									# NM, custom
		sam_line += 'NM:i:%d' % (sam_NM);									# NM, custom
		# sam_line += 'AS:i:%s\t' % (score.strip());			# AS, custom

		return sam_line;

def get_line(fp):
	line = fp.readline();
	if (len(line) == 0):
		return None;
	line = line.replace('\n', '');
	return line;

def complement_base(base):
	if (base == 'A'):
		return 'T';
	if (base == 'C'):
		return 'G';
	if (base == 'T'):
		return 'A';
	if (base == 'G'):
		return 'C';
	return 'N';

def revcomp_seq(sequence):
	ret_seq = '';
	i = 0;
	while (i < len(sequence)):
		ret_seq += complement_base(sequence[len(sequence)-i-1]);
		i += 1;
	return ret_seq;



def convert_to_sam(alignment_file, daligner_reference, daligner_reads, header_conversion_hash, out_sam_file):
	try:
		fp_in = open(alignment_file, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % alignment_file);
		exit(1);

	try:
		fp_out = open(out_sam_file, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % out_sam_file);
		exit(1);

	[ref_headers, ref_lengths, ref_daligner_seq_id] = get_fastq_headers_and_lengths(daligner_reference);
	[read_headers, read_seqs, read_quals] = read_fastq(daligner_reads);
	# [read_headers, read_seqs, read_quals] = read_fastq(daligner_reads);

	STATE_INIT = 0;
	STATE_HEADER = 1;
	STATE_OVERLAP_LINE = 9;
	STATE_ALIGNMENT = 6;
	STATE_MATCHING_LINE = 7;
	STATE_QUERY_LINE = 8;
	STATE_CHECK_IF_NEXT_OVERLAP = 10;

	current_state = STATE_INIT;
	next_state = STATE_INIT;

	# escherichia_coli.fa-dalignerreference.fasta.reads-pacbio-dalignerreads.fasta: 988 records
	p_init = re.compile("\s*%s.%s: ([\d,]+) records\s*" % (os.path.basename(daligner_reference), os.path.basename(daligner_reads)));
	# 1       1 n   [1,356,122..1,358,765] x [   817.. 3,601] :   =    432 diffs  ( 27 trace pts)
	# p_overlap_line = re.compile("(\d+)\s+(\d+)\s+([cn])\s+\[\s*(.*?)\s*\.\.\s*(.*?)\s*\]\s+x\s+\[\s*(.*?)\s*\.\.\s*(.*?)\s*\]\s+:\s+=\s+(\d+)\s+diffs\s*\(\s*(\d+)\s*trace pts\)");
	p_overlap_line = re.compile("\s*([\d,]+)\s+([\d,]+)\s+([cn])\s+\[\s*(.*?)\s*\.\.\s*(.*?)\s*\]\s*x\s*\[\s*(.*?)\s*\.\.\s*(.*?)\s*\]\s*:\s*=\s*([\d,]+)\s*diffs\s*\(\s*([\d,]+)\s*trace pts\)\s*");
	#     1356113 taattaacgc[tgtggcgg-taactaaatcgaagaacagcgccgacaacgcgacaatcccgaccataatga-cgttgag-tgccggagtcc--gccattt
	# p_ref_line = re.compile("(\d+)\s+([actgn\-\[\]\.]+)");
	p_ref_line = re.compile("\s*(\d*)\s*([actgn\-\[\]\. ]+)\s*");
	p_matching_line = re.compile("\s*([:\[\|\]\* ]+)\s*");
	# p_query_line = re.compile("(\d+)\s+([actgn\-\[\]\.]+)\s+(.*?%)");
	p_query_line = re.compile("\s*(\d*)\s*([actgn\-\[\]\.]+)\s*(.*?%){0,1}\s*");

	# print "\s*%s.%s: (\d+) records\s*" % (os.path.basename(daligner_reference), os.path.basename(daligner_reads));

	# Write the SAM header.
	i = 0;
	while i < len(ref_headers):
		line = '@SQ\tSN:%s\tLN:%d\n' % (ref_headers[i], ref_lengths[i]);
		fp_out.write(line);
		i += 1;

	ovl = None;
	num_sam_lines = 0;

	# for line in fp_in:
	while (True):

		if (current_state == STATE_INIT):
			line = get_line(fp_in);
			if (line == None):
				break;

			m = p_init.match(line);
			if (m):
				next_state = STATE_OVERLAP_LINE;
			else:
				next_state = current_state;



		elif (current_state == STATE_OVERLAP_LINE):
			### Empty line.
			line = get_line(fp_in);
			if (line == None):
				break;
			### The actual overlap line.
			line = get_line(fp_in);
			if (line == None):
				break;

 			# if ('1     362 n   [1,075,548..1,079,360]' in line.strip()):
 			# 	print 'Tu sam 2!';

			m = p_overlap_line.match(line);
			if (m):
				ovl = Overlap(line, m);
				next_state = STATE_ALIGNMENT;
	 			# if ('1     362 n   [1,075,548..1,079,360]' in line.strip()):
	 			# 	print 'Tu sam 2.1!';
	 			# 	print ovl.original_line;
			else:
				next_state = current_state;
				sys.stderr.write('ERROR: Daligner line not formatted properly (STATE_OVERLAP_LINE)! Maybe the format changed?\n');
				sys.stderr.write('Line: "%s".\n' % line);
				sys.stderr.write('Exiting.\n');
				exit(1);

		elif (current_state == STATE_ALIGNMENT):
			### Empty line.
			line = get_line(fp_in);
			if (line == None):
				break;

			### Reference line.
			line = get_line(fp_in);
			if (line == None):
				break;
			m = p_ref_line.match(line);

			# if ('1     362 n   [1,075,548..1,079,360]' in ovl.original_line.strip()):
			# 	print 'Tu sam 3!';
			# 	print ovl.verbose_as_string();

			if (m):
				ovl.add_ref_alignment(m.groups()[1]);

			else:
				next_state = current_state;
				sys.stderr.write('ERROR: Daligner line not formatted properly (STATE_ALIGNMENT)! Maybe the format changed?\n');
				sys.stderr.write('Line: "%s".\n' % line);
				sys.stderr.write('Exiting.\n');
				exit(1);

			offset = line.index(m.groups()[1]);
			offset_len = len(m.groups()[1]);

			### Matching line.
			line = get_line(fp_in);
			if (line == None):
				break;
			ovl.add_matching_alignment(line[offset:(offset+offset_len)]);
			# m = p_matching_line.match(line);
			# if (m):
			# 	ovl.add_matching_alignment(m.groups()[0]);
			# else:
			# 	next_state = current_state;
			# 	sys.stderr.write('ERROR: Daligner line not formatted properly (STATE_MATCHING_LINE)! Maybe the format changed?\n');
			# 	sys.stderr.write('Line: "%s".\n' % line);
			# 	sys.stderr.write('Exiting.\n');
			# 	exit(1);

			### Query line.
			line = get_line(fp_in);
			if (line == None):
				break;
			ovl.add_query_alignment(line[offset:(offset+offset_len)]);
			# m = p_query_line.match(line);
			# if (m):
			# 	ovl.add_query_alignment(m.groups()[1]);
			# else:
			# 	next_state = current_state;
			# 	sys.stderr.write('ERROR: Daligner line not formatted properly (STATE_QUERY_LINE)! Maybe the format changed?\n');
			# 	sys.stderr.write('Line: "%s".\n' % line);
			# 	sys.stderr.write('Exiting.\n');
			# 	exit(1);

			next_state = STATE_CHECK_IF_NEXT_OVERLAP;

		elif (current_state == STATE_CHECK_IF_NEXT_OVERLAP):
			fp_tell = fp_in.tell();
			### Empty line.
			line = get_line(fp_in);
			if (line == None):
				break;

			### Either a new overlap or a continuation to the previous alignment.
			line = get_line(fp_in);
			if (line == None):
				break;

			### Return to the previous position, because it needs to be re-read.
			fp_in.seek(fp_tell, 0);

 # 1     362 n   [1,075,548..1,079,360] x [     0.. 4,045] :   =    670 diffs  ( 39 trace pts)
 			# if ('1     362 n   [1,075,548..1,079,360]' in ovl.original_line.strip()):
 			# 	print 'Tu sam 1!';

			m = p_overlap_line.match(line);
			if (m):
				sam_line = ovl.convert_to_sam(ref_headers, ref_daligner_seq_id, read_headers, read_seqs, read_quals, header_conversion_hash);
				num_sam_lines += 1;
				fp_out.write(sam_line + '\n');
	 			# if ('1     362 n   [1,075,548..1,079,360]' in ovl.original_line.strip()):
					# print 'Tu sam 4!';
					# print ovl.verbose_as_string();

				next_state = STATE_OVERLAP_LINE;
			else:
				next_state = STATE_ALIGNMENT;

		current_state = next_state;

	if (ovl):
		sam_line = ovl.convert_to_sam(ref_headers, ref_daligner_seq_id, read_headers, read_seqs, read_quals, header_conversion_hash);
		num_sam_lines += 1;
		fp_out.write(sam_line + '\n');

	fp_in.close();
	fp_out.close();

	sys.stderr.write('Number of outputted SAM lines: %d\n' % (num_sam_lines));



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#	reads_file			Path to a FASTA/FASTQ file containing reads.
#	reference_file		Path to a reference genome FASTA file.
#	machine_name		A symbolic name to specify a set of parameters for a specific sequencing platform.
#	output_path			Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#	output_suffix		A custom suffix that can be added to the output filename.
def run(run_type, reads_file, reference_file, machine_name, output_path, output_suffix=''):
	parameters = '';
	num_threads = multiprocessing.cpu_count() / 2;

	if ((machine_name.lower() == 'illumina') or (machine_name.lower() == 'roche')):
		# parameters = '-v -s1 -h10 -e.9';
		### I get poor results on Illumina data (simulated), concretely DALIGNER mapps 0 reads. I think the problem is 'alignment but
		### simply a set of trace points, typically every 100bp or so, that allow the', and reads that I simulated were 150bp in length.
		parameters = '-v';ŠP

	elif ((machine_name.lower() == 'pacbio')):
		# parameters = '-t %s -x pacbio' % str(num_threads);
		parameters = '-v';

	elif ((machine_name.lower() == 'nanopore')):
		parameters = '-v -e.7 -k10';

	elif ((machine_name.lower() == 'k9')):
		parameters = '-v -e.7 -k9';

	elif ((machine_name.lower() == 'k10')):
		parameters = '-v -e.7 -k10';

	# elif ((machine_name.lower() == 'debug')):
	# 	parameters = '-t %s' % str(num_threads);

	else:			# default
		parameters = '-vd';



	if (output_suffix != ''):
		if (output_suffix.lower().endswith('.sam')):
			output_filename = os.path.splitext(output_suffix)[0];
		else:
			output_filename = '%s-%s' % (MAPPER_NAME, output_suffix);
	else:
		output_filename = MAPPER_NAME;
	
	# Check if the given input file is a FASTA or FASTQ, and convert to FASTA if necessary.
	if (reads_file[-1] == 'q'):
		sys.stderr.write('[%s wrapper] Converting FASTQ to FASTA...\n' % (MAPPER_NAME));
		reads_fasta = reads_file[0:-1] + 'a';
		convert_to_fasta(reads_file, reads_fasta);
		reads_file = reads_fasta;
		sys.stderr.write('\n');

	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	sam_file = '%s/%s.sam' % (output_path, output_filename);
	memtime_file = '%s/%s.memtime' % (output_path, output_filename);
	memtime_file_hpcmapper = '%s/%s-hpcmapper.memtime' % (output_path, output_filename);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, output_filename);

	### Convert the input files to absolute paths.
	if (os.path.isabs(reads_file) == False):
		reads_file = os.path.abspath(reads_file);
	if (os.path.isabs(reference_file) == False):
		reference_file = os.path.abspath(reference_file);

	daligner_reference_file = reference_file + '-dalignerreference.fasta';
	if (run_type == 'align' or run_type == 'run'):
		index_file = daligner_reference_file + '.dam';

		if (os.path.exists(index_file)):
			sys.stderr.write('[%s wrapper] The DALIGNER index file already exists ("%s"), removing.\n' % (MAPPER_NAME, index_file));
			os.remove(index_file);

		# Run the indexing process, and measure execution time and memory.
		# daligner_reference_file = reference_file if (reference_file.lower().endswith('fasta')) else (reference_file + '.fasta');
		sys.stderr.write('[%s wrapper] Wrapping the sequences in the reference FASTA file. DALIGNER has a line length limit of 9998 chars.\n' % (MAPPER_NAME));
		wrap_fasta_file(reference_file, daligner_reference_file);

		if (True or (not os.path.exists(index_file))):
			if (not os.path.exists(daligner_reference_file)):
				sys.stderr.write('[%s wrapper] Copying reference to satisfy the extension requirements...\n' % (MAPPER_NAME));
				command = 'cp %s %s.fasta' % (reference_file, reference_file);
				execute_command(command);

			sys.stderr.write('[%s wrapper] Generating index...\n' % (MAPPER_NAME));
			command = '%s %s/fasta2DAM %s %s' % (measure_command_wrapper(memtime_file_index), ALIGNER_DB_PATH, index_file, daligner_reference_file);
			execute_command(command);
			sys.stderr.write('\n');
		else:
			sys.stderr.write('[%s wrapper] Reference index already exists. Continuing.\n' % (MAPPER_NAME));
			sys.stderr.flush();

	daligner_reads_file = '%s-dalignerreads.fasta' % (os.path.splitext(reads_file)[0]);
	if (True or (not os.path.exists(daligner_reads_file))):
		sys.stderr.write('[%s wrapper] Modifying the reads file to have PacBio headers...\n' % (MAPPER_NAME));
		# command = 'cp %s %s.fasta' % (reads_file, reads_file);
		# subprocess.call(command, shell=True);
		header_conversion_hash = convert_reads_to_pacbio_format(reads_file, daligner_reads_file);
		sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] Converting the reads file into a DB file...\n' % (MAPPER_NAME));
	daligner_reads_file_db = '%s.db' % (daligner_reads_file);
	if (os.path.exists(daligner_reads_file_db)):
		sys.stderr.write('[%s wrapper] The DALIGNER reads DB file already exists ("%s"), removing.\n' % (MAPPER_NAME, daligner_reads_file_db));
		os.remove(daligner_reads_file_db);
	command = '%s %s/fasta2DB %s.db %s' % (measure_command_wrapper(memtime_file_index), ALIGNER_DB_PATH, daligner_reads_file, daligner_reads_file);
	execute_command(command);
	sys.stderr.write('\n');

	### DALIGNER's HPCmapper script basically just generates a shell script with commands that need to be run to generate alignments in parallel and then join them into one LAS file.
	### Instead of outputting this script to a file, we intercept the STDOUT and modify it a bit.
	### Modifications are needed because DALIGNER's generated script expects that it's binaries are in PATH, and also it generates intermediate files in the current folder.
	### That's why we modify the PATH variable first, and change the execution folder to the output folder.
	### Please note that reads_file and reference_file then need to be absolute paths, so that's why we performed the conversion above.
	# Run the alignment process, and measure execution time and memory.	
	#
	#
	#
	### Extracting alignments is not documented very well, or at least there is not any examples out there on how to do that.
	### From the DALIGNER's README on GitHub, it says:
	### "If the -c option is given then a cartoon rendering is displayed, and if -a or -r option is set then an alignment of the local alignment is displayed."
	###
	### (1) The -c on its own doesn't really show much useful info, e.g.:
	###	  1       1 n   [1,356,122..1,358,765] x [   817.. 3,601]  ( 27 trace pts)
	###        1356122                3280910
	###    A ==========+------------+=========>          dif/(len1+len2) = 432/(2643+2784) = 15.92%
	###    B     ======+------------>
	###            817 
	###
	### (2) The -a option produces alignments in the BLAST like format, that is hard to parse. E.g.:
	### 1       1 n   [1,356,122..1,358,765] x [   817.. 3,601] :   =    432 diffs  ( 27 trace pts)
	###
	###     1356113 taattaacgc[tgtggcgg-taactaaatcgaagaacagcgccgacaacgcgacaatcccgaccataatga-cgttgag-tgccggagtcc--gccattt
	###             ::::::::::[||||||||*|||||||||||*|| *|||||||||||||||*|*|||||||||||||*||||*|||||||*|||||||||||**|||||||
	###         808 tagttcacct[tgtggcgggtaactaaatcgtag-acagcgccgacaacgtg-caatcccgaccattatgaacgttgaggtgccggagtccctgccattt  11.2%
	###
	### (3) The -r option - causes slightly different row width.
	### 1       1 n   [1,356,122..1,358,765] x [   817.. 3,601] :   =    432 diffs  ( 27 trace pts)
	###
	###     1356113 taattaacgc[tgtggcgg-taactaaatcgaagaacagcgccgacaacgcgacaatcccgaccataatga-cgttgag-tgccggagtcc--g
	###             ::::::::::[||||||||*|||||||||||*||*|||||||||||||||*|*|||||||||||||*||||*|||||||*|||||||||||**|
	###         808 tagttcacct[tgtggcgggtaactaaatcgtag-acagcgccgacaacgtg-caatcccgaccattatgaacgttgaggtgccggagtccctg  12.0%
	sys.stderr.write('[%s wrapper] Running %s...\n' % (MAPPER_NAME, MAPPER_NAME));
	if (run_type == 'align' or run_type == 'run'):
		command = '%s %s/HPCmapper %s %s %s' % (measure_command_wrapper(memtime_file_hpcmapper), ALIGNER_PATH, parameters, daligner_reference_file, daligner_reads_file);
		[out, err] = execute_command_get_stdout(command);

		print out;

		### Replace ampersands with '\n' so it's easier to split commands and add measurement calls to the commands.
		out = out.replace('&&', '\n');
		### LAshow should extract the overlaps/alignments from the LAS file.
		las_file = '%s.%s.las' % (os.path.basename(daligner_reference_file), os.path.basename(daligner_reads_file));
		out += '\nLAshow -a %s %s %s > %s.txt' % (daligner_reference_file, daligner_reads_file, las_file, las_file);
		# commands_daligner = 'PATH="$PATH:%s"\necho $PATH\ncd %s\n%s %s\nLAshow -a %s %s %s > %s.txt' % (ALIGNER_PATH, output_path, measure_command_wrapper(memtime_file), out, daligner_reference_file, daligner_reads_file, las_file, las_file);

		# for line in out.split('\n'):
		# 	print '%s\n' % line;
		# print '---------------';

		### Prepare measurement command for each line of the generated script.
		daligner_out_formatted = [command.split('#')[0].strip() for command in out.split('\n') if (len(command.split('#')[0].strip()) > 0)]
		memtime_files = [];
		i = 0;
		while (i < len(daligner_out_formatted)):
			memtime_file_i = '%s/%s-%d.memtime' % (output_path, output_filename, i);
			memtime_files.append(memtime_file_i);
			measure_command_i = measure_command_wrapper(memtime_file_i);
			daligner_out_formatted[i] = '%s %s' % (measure_command_i, daligner_out_formatted[i]);
			i += 1;
		joined_commands = '; '.join(daligner_out_formatted);

		### Add the paths so the commands can be executed properly.
		commands_daligner = 'PATH="$PATH:%s"; echo $PATH; cd %s; %s' % (ALIGNER_PATH, output_path, joined_commands);
		execute_command(commands_daligner);
		sys.stderr.write('\n');

		final_memtime_file = memtime_file;
		parse_memtime_files_and_accumulate(memtime_files, final_memtime_file);

		for memtime_file in memtime_files:
			try:
				os.remove(memtime_file);
			except OSError:
				pass;

		sys.stderr.write('[%s wrapper] Converting the output to SAM format at path "%s"...\n' % (MAPPER_NAME, sam_file));
		convert_to_sam('%s/%s.txt' % (output_path, las_file), daligner_reference_file, daligner_reads_file, header_conversion_hash, sam_file);

	elif (run_type == 'overlap'):
		command = '%s %s/HPCdaligner %s %s' % (measure_command_wrapper(memtime_file_hpcmapper), ALIGNER_PATH, parameters, daligner_reads_file);
		[out, err] = execute_command_get_stdout(command);
		### LAshow should extract the overlaps/alignments from the LAS file.
		las_file = '%s.1.las' % (daligner_reads_file);
		commands_daligner = 'PATH="$PATH:%s"\necho $PATH\ncd %s\n%s\nLAshow -a %s %s %s > %s.txt' % (ALIGNER_PATH, output_path, out, daligner_reference_file, daligner_reads_file, las_file, las_file);
		commands_daligner = '; '.join([command for command in commands_daligner.split('\n') if (len(command) > 0 and command[0] != '#')]);
		execute_command(commands_daligner);
		sys.stderr.write('\n');

	elif (run_type == 'onlyconvert'):
		sys.stderr.write('[%s wrapper] Converting the output to SAM format at path "%s"...\n' % (MAPPER_NAME, sam_file));
		las_file = '%s.%s.las' % (os.path.basename(daligner_reference_file), os.path.basename(daligner_reads_file));
		convert_to_sam('%s/%s.txt' % (output_path, las_file), daligner_reference_file, daligner_reads_file, header_conversion_hash, sam_file);

	sys.stderr.write('[%s wrapper] %s wrapper script finished processing.\n' % (MAPPER_NAME, MAPPER_NAME));

	return sam_file


# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
	sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (MAPPER_NAME, MAPPER_NAME));
	sys.stderr.write('[%s wrapper] Creating a folder for all %s repos...\n' % (MAPPER_NAME, MAPPER_NAME));
	command = 'mkdir -p %s/%s' % (ALIGNERS_PATH_ROOT_ABS_, MAPPER_NAME);
	execute_command(command);

	sys.stderr.write('[%s wrapper] Cloning git repository.\n' % (MAPPER_NAME));
	command = 'cd %s/%s; git clone %s' % (ALIGNERS_PATH_ROOT_ABS_, MAPPER_NAME, ALIGNER_URL);
	execute_command(command);

	sys.stderr.write('[%s wrapper] Running make.\n' % (MAPPER_NAME));
	command = 'cd %s; make' % (ALIGNER_PATH);
	execute_command(command);

	sys.stderr.write('[%s wrapper] Cloning git repository.\n' % (MAPPER_NAME));
	command = 'cd %s/%s; git clone %s' % (ALIGNERS_PATH_ROOT_ABS_, MAPPER_NAME, ALIGNER_DB_URL);
	execute_command(command);

	sys.stderr.write('[%s wrapper] Running make.\n' % (MAPPER_NAME));
	command = 'cd %s; make' % (ALIGNER_DB_PATH);
	execute_command(command);

	# sys.stderr.write('[%s wrapper] Checking out commit "eb428d7d31ced059ad39af2701a22ebe6d175657" for reproducibility purposes.\n' % (MAPPER_NAME));
	# command = 'cd %s; git checkout eb428d7d31ced059ad39af2701a22ebe6d175657' % (ALIGNER_PATH);
	# subprocess.call(command, shell='True');
	# sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] All installation steps finished.\n' % (MAPPER_NAME));
	sys.stderr.write('\n');



def verbose_usage_and_exit():
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode [<reads_file> <reference_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\t- mode          - either "align", "overlap" or "install". If "install" other parameters can be ommitted.\n');
	sys.stderr.write('\t- machine_name  - "illumina", "roche", "pacbio", "nanopore" or "default".\n');
	sys.stderr.write('\t- output_suffix - suffix for the output filename. If this parameter ends with ".sam", the value will be used as full output filename.\n');

	exit(0);

if __name__ == "__main__":
	if (len(sys.argv) < 2 or len(sys.argv) > 7):
		verbose_usage_and_exit();

	if (sys.argv[1] == 'install'):
		download_and_install();
		exit(0);

	elif (sys.argv[1] == 'align' or sys.argv[1] == 'run' or sys.argv[1] == 'onlyconvert'):
		if (len(sys.argv) < 6):
			verbose_usage_and_exit();

		reads_file = sys.argv[2];
		reference_file = sys.argv[3];
		machine_name = sys.argv[4];
		output_path = sys.argv[5];
		output_suffix = '';

		if (len(sys.argv) == 7):
			output_suffix = sys.argv[6];
		run(sys.argv[1], reads_file, reference_file, machine_name, output_path, output_suffix);

	elif (sys.argv[1] == 'overlap'):
		if (len(sys.argv) < 6):
			verbose_usage_and_exit();

		reads_file = sys.argv[2];
		reference_file = sys.argv[3];
		machine_name = sys.argv[4];
		output_path = sys.argv[5];
		output_suffix = '';

		if (len(sys.argv) == 7):
			output_suffix = sys.argv[6];
		run(sys.argv[1], reads_file, reference_file, machine_name, output_path, output_suffix);

	else:
		verbose_usage_and_exit();
