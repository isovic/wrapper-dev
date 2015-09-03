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

# ALIGNER_URL = 'https://github.com/isovic/marginAlign.git'; # 'git://github.com/benedictpaten/marginAlign.git';
ALIGNER_URL = 'git://github.com/benedictpaten/marginAlign.git';
ALIGNER_PATH = ALIGNERS_PATH_ROOT_ABS_ + '/marginAlign/';
BIN = 'marginAlign';
MAPPER_NAME = 'marginAlign';

RUNNING_PATH = os.path.dirname(sys.argv[0]);



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

def read_fastq(fastq_path):
	headers = [];
	seqs = [];
	quals = [];
	
	fp_in = None;

	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % fastq_path);
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

def check_if_fastq(fastq_path):
	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % fastq_path);
		exit(1);
	[header, read] = get_single_read(fp_in);
	fp_in.close();
	if (len(read) == 4):
		return True;
	return False;



# def wrap_fasta_file(fasta_file, daligner_fasta_file):
# 	try:
# 		fp_in = open(fasta_file, 'r');
# 	except:
# 		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % fasta_file);
# 		exit(0);

# 	try:
# 		fp_out = open(daligner_fasta_file, 'w');
# 	except:
# 		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % daligner_fasta_file);
# 		exit(0);

# 	current_read = 0;

# 	while True:
# 		[header, read] = get_single_read(fp_in);
		
# 		if (len(read) == 0):
# 			break;

# 		current_read += 1;

# 		if (len(read[1]) <= 10):	### DALIGNER has a lower length limit of 10bp.
# 			continue;

# 		read[1] = re.sub("(.{500})", "\\1\n", read[1], 0, re.DOTALL);	### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
# 		if (len(read) == 4):
# 			read[3] = re.sub("(.{500})", "\\1\n", read[3], 0, re.DOTALL);	### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
# 		fp_out.write('\n'.join(read) + '\n');

# 	sys.stderr.write('\n');
# 	fp_in.close();
# 	fp_out.close();

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


### Modifies the FASTQ headers to replace all special chars with '_'.
def modify_reference_headers(input_fastq_path, out_fastq_path):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);
	try:
		fp_out = open(out_fastq_path, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % out_fastq_path);
		exit(0);
	num_matches = 0;
	header_hash = {};
	while True:
		[header, read] = get_single_read(fp_in);
		if (len(read) == 0):
			break;
		new_header = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
		header_hash[new_header[1:]] = read[0][1:];
		read[0] = new_header;
		fp_out.write('\n'.join(read) + '\n');
	sys.stderr.write('\n');
	fp_in.close();
	fp_out.close();
	return header_hash;

### Modifies the FASTQ headers to replace all special chars with '_'.
### Also, filters out sequences longer than 50000bp, because marginAlign passes CIGAR strings as space-separated ops through arguments to another program.
### Since there are limitations to the number of arguments to a program given through shell (at least when Python's subprocess package is used to address the shell),
### the program crashes.
def modify_read_headers_and_remove_long_ones(input_fastq_path, out_fastq_path):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);
	try:
		fp_out = open(out_fastq_path, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % out_fastq_path);
		exit(0);
	num_matches = 0;
	header_hash = {};
	while True:
		[header, read] = get_single_read(fp_in);
		if (len(read) == 0):
			break;
		if (len(read[1]) <= 50000):
			# read[0] = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
			new_header = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
			header_hash[new_header[1:]] = read[0][1:];
			read[0] = new_header;
			fp_out.write('\n'.join(read) + '\n');
	sys.stderr.write('\n');
	fp_in.close();
	fp_out.close();
	return header_hash;

def fix_sam_qnames_after_marginAlign(input_sam_path, ref_header_hash, read_header_hash, out_sam_path):
	if (input_sam_path == out_sam_path):
		sys.stderr.write('ERROR: Input and output SAM files are the same! Skipping the update of qname and rname values.\n');
		return;
	try:
		fp_in = open(input_sam_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_sam_path);
		exit(0);
	try:
		fp_out = open(out_sam_path, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % out_sam_path);
		exit(0);
	i = 0;
	for line in fp_in:
		i += 1;
		line = line.strip();
		if (len(line) == 0 or line[0] == '@'):
			if (line.startswith('@SQ')):
				split_line = line.split('\t');
				found_hname = False;
				for param in split_line:
					if (param.startswith('SN:')):
						hname = param.split(':')[-1];
						try:
							original_hname = ref_header_hash[hname];
						except:
							original_hname = hname;
						new_line = line.replace(hname, original_hname);
						fp_out.write(new_line + '\n');
						found_hname = True;
						break;
				if (found_hname == False):
					fp_out.write(line + '\n');

			else:
				fp_out.write(line + '\n');
			continue;
		sys.stderr.write('\rLine %d' % (i));

		split_line = line.split('\t');
		qname = split_line[0];
		rname = split_line[2];

		if (qname == '*' or rname == '*'):
			fp_out.write(line + '\n');
			continue;

		try:
			original_qname = read_header_hash[qname];
		except:
			original_qname = qname;
		try:
			original_rname = ref_header_hash[rname];
		except:
			original_rname = qname;

		new_line = line.replace(qname, original_qname);
		new_line = new_line.replace(rname, original_rname);
		fp_out.write(new_line + '\n');
	sys.stderr.write('\n');

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
		parameters = ' ';

	elif ((machine_name.lower() == 'pacbio')):
		# parameters = '-t %s -x pacbio' % str(num_threads);
		parameters = ' ';

	elif ((machine_name.lower() == 'nanopore')):
		parameters = ' ';

	elif ((machine_name.lower() == 'nanopore1d')):
		parameters = ' ';

	elif ((machine_name.lower() == 'nanopore2d')):
		parameters = ' ';

	# elif ((machine_name.lower() == 'debug')):
	# 	parameters = '-t %s' % str(num_threads);

	else:			# default
		parameters = ' --em';

	used_mapper = 'last';

	if (output_suffix != ''):
		if (output_suffix.lower().endswith('.sam')):
			output_filename = os.path.splitext(output_suffix)[0];
		else:
			output_filename = '%s-%s' % (MAPPER_NAME, output_suffix);
	else:
		output_filename = MAPPER_NAME;
	
	### Check if the given input file is a FASTA or FASTQ, and convert to FASTA if necessary.
	if (check_if_fastq(reads_file) == False):
		sys.stderr.write('[%s wrapper] marginAlign requires input reads to be in the FASTQ format. Qualities are missing! Exiting.\n' % (MAPPER_NAME));
		exit(1);

	### Convert the input files to absolute paths.
	if (os.path.isabs(reads_file) == False):
		reads_file = os.path.abspath(reads_file);
	if (os.path.isabs(reference_file) == False):
		reference_file = os.path.abspath(reference_file);

	### marginAlign has some quirks with input files. For example, the headers cannot contain special characters (such as '|') because
	### marginAlign crashes. Since NCBI headers are heavy on those, we need to make a temporary reads file which will modify the headers.
	sys.stderr.write('[%s wrapper] Creating a copy of the reference file with fixed headers.\n' % (MAPPER_NAME));
	marginAlign_reference_file = os.path.splitext(reference_file)[0] + '-marginAlign.fa';
	ref_header_hash = modify_reference_headers(reference_file, marginAlign_reference_file);

	### The special chars can also make a problem if present in the reads headers. But additionally, if the reads are too long,
	### marginAlign will crash because it tries to pass the SAM alignments to the PecanRealign via commandline arguments.
	sys.stderr.write('[%s wrapper] Creating a copy of the reads file with fixed headers and maximum read length of 50000bp.\n' % (MAPPER_NAME));
	marginAlign_reads_file = os.path.splitext(reads_file)[0] + '-marginAlign.fastq';
	read_header_hash = modify_read_headers_and_remove_long_ones(reads_file, marginAlign_reads_file);

	temp_sam_file = '%s/%s.sam' % (output_path, output_filename);
	sam_file = '%s/fixed-%s.sam' % (output_path, output_filename);
	memtime_file = '%s/%s.memtime' % (output_path, output_filename);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, output_filename);

	# output_model_file = '%s/hmm-%s-%s.txt' % (reads_path, reads_basename, used_mapper);
	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	reads_path = os.path.dirname(reads_file);
	reference_path = os.path.dirname(reference_file);
	output_model_file = '%s/hmm-%s-%s.txt' % (reference_path, machine_name.lower(), used_mapper);
	jobtree = '%s/jobTree' % (output_path);
	

	if (run_type == 'align' or run_type == 'run'):
		sys.stderr.write('[%s wrapper] Running %s...\n' % (MAPPER_NAME, MAPPER_NAME));
		if ((not os.path.exists(output_model_file)) or ('--em' in parameters)):
			sys.stderr.write('[%s wrapper] Expectation maximization will be used.\n' % (MAPPER_NAME));
			if (('--em' in parameters) == False):
				parameters += ' --em';
			if (('--outputModel' in parameters) == False):
				parameters += ' --outputModel=%s' % output_model_file;
		else:
			sys.stderr.write('[%s wrapper] Existing model found. Will be using the model values.\n' % (MAPPER_NAME));
			parameters += ' --inputModel=%s' % output_model_file;

		if (os.path.exists(jobtree)):
			sys.stderr.write('[%s wrapper] Removing old jobtree folder.\n' % (MAPPER_NAME));
			execute_command('rm -r %s' % (jobtree));

		# execute_command('%s %s/marginAlign %s %s %s --jobTree %s %s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), ALIGNER_PATH, marginAlign_reads_file, marginAlign_reference_file, temp_sam_file, jobtree, parameters, num_threads, num_threads));

	sys.stderr.write('[%s wrapper] Fixing SAM qname and rname headers to original values. Temp SAM file: "%s", final SAM file: "%s".\n' % (MAPPER_NAME, temp_sam_file, sam_file));
	fix_sam_qnames_after_marginAlign(temp_sam_file, ref_header_hash, read_header_hash, sam_file);

	sys.stderr.write('[%s wrapper] %s wrapper script finished processing.\n' % (MAPPER_NAME, MAPPER_NAME));

	return sam_file


# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
	sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (MAPPER_NAME, MAPPER_NAME));

	sys.stderr.write('[%s wrapper] Cloning git repository.\n' % (MAPPER_NAME));
	command = 'cd %s; git clone %s' % (ALIGNERS_PATH_ROOT_ABS_, ALIGNER_URL);
	execute_command(command);
# 
	sys.stderr.write('[%s wrapper] Initializing submodules.\n' % (MAPPER_NAME));
	command = 'cd %s; git submodule update --init' % (ALIGNER_PATH);
	execute_command(command);

	sys.stderr.write('[%s wrapper] Running make.\n' % (MAPPER_NAME));
	command = 'cd %s; make' % (ALIGNER_PATH);
	execute_command(command);

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

	elif (sys.argv[1] == 'align' or sys.argv[1] == 'run'):
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
