#!/usr/local/bin/python3

# Michael Matschiner, 2015-03-26
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(1)
import argparse, textwrap, random, tempfile, os, subprocess, statistics, re

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      This script uses one fasta format query file with one or
      more sequence entries and one or more subject sequence
      files that must also be in fasta format, and must have
      blast databases prepared for them with makeblastdb.
      Putative orthologs to query sequences are detected in
      subject sequences based on variable settings for
      strictness and evalue and/or bitscore thresholds. The
      query sequence(s) can be divided into sliding windows.
      Blast results are used to build a (multiple) sequence
      alignment and exported in fasta or phylip format.
      Note that blast+ must be installed on your system for
      this script to run (also, blast+ 2.2.26 and 2.2.30 work,
      but 2.2.29 throws an error when the query is composed of
      only missing data). On the other hand, blast+ 2.2.26
      is not compatible with genetic codes other than the
      standard code.
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.77',
    help='show program version number and exit.'
    )
parser.add_argument(
    'query',
    nargs=1,
    type=argparse.FileType('r'),
    help='the query file name (fasta).'
    )
parser.add_argument(
    'subject',
    nargs='+',
    type=str,
    help="one or more blast database names. If a single file name is provided \
    and the first character of the file is not '>', it is assumed to be a list \
    of blast database names.")
parser.add_argument(
    '-s', '--strictness',
    nargs=1,
    type=str,
    default=[0],
    metavar="INT",
    dest='strictness',
    help="strictness setting for ortholog identification. \
    0: use all subjects with sufficient fit.\
    1: per query fragment region, use the best subject with sufficient fit.\
    2: as 1, but all subjects must be on the same scaffold/contig, have the \
    same orientation, and be ordered the same way in the query and the subject."
    )
parser.add_argument(
    '-e', '--evalue',
    nargs=1,
    type=float,
    default=[None],
    metavar="FLOAT",
    dest='evalue',
    help="evalue threshold to consider blast hits (regardless of specified \
    evalue, only the best blast hit per alignment region is used). Is \
    overridden for query sequences that include '[&evalue=FLOAT]' (where \
    FLOAT is a floating point number) in their sequence id."
    )
parser.add_argument(
    '-b', '--bitscore',
    nargs=1,
    type=float,
    default=[None],
    metavar="FLOAT",
    dest='bitscore',
    help="bitscore threshold to consider blast hits (regardless of specified \
    bitscore, only the best blast hit per alignment region is used). Is \
    overridden for query sequences that include '[&bitscore=FLOAT]' (where \
    FLOAT is a floating point number) in their sequence id."
    )
parser.add_argument(
    '-t', '--translate',
    action='store_true',
    dest='translate',
    help="translate subject sequences (use with amino acid sequence query)."
    )
parser.add_argument(
    '-c', '--genetic-code',
    nargs=1,
    type=float,
    default=[None],
    metavar="INT",
    dest='genetic_code',
    help="If amino acid sequence queries are used as queries for translated \
    subject sequences, the genetic code can be specified with this option. \
    Note that at least blast+ 2.2.29 is required with this option.\
    See here for a list of all available genetic code:\
    http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?"
    )
parser.add_argument(
    '-wl', '--window-length',
    nargs=1,
    type=int,
    default=[None],
    metavar="INT",
    dest='window_length',
    help="length of the sliding window."
    )
parser.add_argument(
    '-ws', '--window-shift',
    nargs=1,
    type=int,
    default=[None],
    metavar="INT",
    dest='window_shift',
    help="sliding window shift distance in bp."
    )
parser.add_argument(
    '-r', '--refine',
    action='store_true',
    dest='refine',
    help="refine the alignment. Note that the software mafft must be installed\
    in order to use this option."
    )
parser.add_argument(
    '-m', '--minimum-completeness',
    nargs=1,
    type=int,
    default=[None],
    metavar="INT",
    dest='minimum_completeness',
    help="minimum number of alignment sites without missing data."
    )
parser.add_argument(
    '--write-empty-files',
    action='store_true',
    dest='write_empty_files',
    help="write empty files if alignment is not sufficiently complete.\
    By default, no files are written for these alignments."
    )
parser.add_argument(
    '-of', '--output-format',
    nargs=1,
    type=str,
    default=["fasta"],
    metavar="STR",
    dest='output_format',
    help="output alignment format: fasta (default) or phylip."
    )
parser.add_argument(
    '--overwrite',
    action='store_true',
    dest='overwrite',
    help="overwrite existing alignments."
    )
args = parser.parse_args()
query = args.query[0]

if len(args.subject) == 1:
    # Check if the file name specified for the subject is a list of
    # subject file names.
    with open(args.subject[0], "r") as f:
        first_line = f.readline()
    if first_line[0:1] == ">":
        subjects = args.subject
    else:
        subject_list_file = open(args.subject[0], "r")
        subject_list = subject_list_file.readlines()
        subjects = []
        for subject_line in subject_list:
            subjects.append(subject_line.strip())
        subject_list_file.close()
else:
    subjects = args.subject
# Make sure that all subjects files are present.
for subject in subjects:
    if not os.path.isfile(subject):
        print("ERROR: Subject file " + subject + " could not be found!")
        sys.exit(1)
    elif not os.path.isfile(subject + ".nhr") and not os.path.isfile(subject + ".00.nhr"):
        print("ERROR: File " + subject + ".nhr for subject " + subject + " could not be found!")
        sys.exit(1)
    elif not os.path.isfile(subject + ".nin") and not os.path.isfile(subject + ".00.nin"):
        print("ERROR: File " + subject + ".nin for subject " + subject + " could not be found!")
        sys.exit(1)
    elif not os.path.isfile(subject + ".nsq") and not os.path.isfile(subject + ".00.nsq"):
        print("ERROR: File " + subject + ".nsq for subject " + subject + " could not be found!")
        sys.exit(1)

strictness = int(args.strictness[0])
if args.evalue[0] != None:
    evalue = float(args.evalue[0])
else:
    evalue = None
if args.bitscore[0] != None:
    bitscore = float(args.bitscore[0])
else:
    bitscore = None
if args.genetic_code[0] != None:
    genetic_code = int(args.genetic_code[0])
else:
    genetic_code = None
translate = args.translate
window_length = args.window_length[0]
window_shift = args.window_shift[0]
if window_shift == None and window_length != None:
    window_shift = window_length
refine = args.refine
minimum_completeness = args.minimum_completeness[0]
write_empty_files = args.write_empty_files
output_format = args.output_format[0]
overwrite = args.overwrite

# Read the query fasta string and store all entries in memory.
sys.stdout.flush()
query_string = query.read()
query_entries = query_string.split(">")
query_ids = []
query_seqs = []
query_evalues = []
query_bitscores = []
pattern = re.compile('\[&(.+)\]')
for query_entry in query_entries:
    if query_entry != "":
        query_entry_lines = query_entry.split("\n")
        pattern_match = pattern.search(query_entry_lines[0])
        if pattern_match == None:
            query_ids.append(query_entry_lines[0])
            query_evalues.append(evalue)
            query_bitscores.append(bitscore)
        else:
            query_ids.append(query_entry_lines[0].replace(pattern_match.group(0),""))
            id_command = pattern_match.group(1)
            evalue_pattern = re.compile('evalue=([\d\.]+)')
            evalue_match = evalue_pattern.search(id_command)
            if evalue_match == None:
                query_evalues.append(evalue)
            else:
                query_evalues.append(float(evalue_match.group(1)))
            bitscore_pattern = re.compile('bitscore=([\d\.]+)')
            bitscore_match = bitscore_pattern.search(id_command)
            if bitscore_match == None:
                query_bitscores.append(bitscore)
            else:
                query_bitscores.append(float(bitscore_match.group(1)))
        query_seqs.append("".join(query_entry_lines[1:]))

# Make sure the query sequences do not contain terminal stop codons. Terminal stop codons are accepted but removed.
for x in range(len(query_seqs)):
    if query_seqs[x][-1] == "*":
        query_seqs[x] = query_seqs[x][:-2]
    if "*" in query_seqs[x]:
        print("ERROR: Non-terminal stop codon found in query sequence "+ query_ids[x] + "! These should be removed.")
        sys.exit(1)

# Feedback.
print("")
print("find_orthologs.py")
print("")
print("Settings:")
if len(query_ids) == 1:
    print("    query:                " + query.name + " (" + str(len(query_ids)) + " sequence)")
else:
    print("    query:                " + query.name + " (" + str(len(query_ids)) + " sequences)")
for subject in subjects:
    print("    subject:              " + subject)
if strictness == 0:
    print("    strictness:           0")
elif strictness == 1:
    print("    strictness:           1")
elif strictness == 2:
    print("    strictness:           2")
else:
    print("ERROR: strictness must be 0, 1, or 2!")
    sys.exit(1)
if evalue == None:
    print("    evalue:               -")
else:
    print("    evalue:               " + str(evalue))
if bitscore == None:
    print("    bitscore:             -")
else:
    print("    bitscore:             " + str(bitscore))
if translate == True:
    print("    translate:            yes")
else:
    print("    translate:            no")
if translate == False:
    print("    genetic_code:         -")
else:
    if genetic_code == None:
        print("    genetic_code:         standard")
    elif genetic_code == 1:
        print("    genetic_code:         standard")
    elif genetic_code == 2:
        print("    genetic_code:         vertebrate mitochondrial")
    elif genetic_code == 3:
        print("    genetic_code:         yeast mitochondrial")
    elif genetic_code == 5:
        print("    genetic_code:         invertebrate mitochondrial")
    else:
        print("    genetic_code:         other")
if window_length == None:
    print("    window_length:        -")
else:
    print("    window_length:        " + str(window_length))
if window_shift == None:
    print("    window_shift:         -")
else:
    print("    window_shift:         " + str(window_shift))
if refine:
    print("    refine:               yes")
else:
    print("    refine:               no")
if minimum_completeness == None:
    print("    minimum_completeness: -")
else:
    print("    minimum_completeness: " + str(minimum_completeness))
if write_empty_files:
    print("    write_empty_files:    yes")
else:
    print("    write_empty_files:    no")
print("    output_format:        " + str(output_format))
if overwrite:
    print("    overwrite:            yes")
else:
    print("    overwrite:            no")
print("")
sys.stdout.flush()

# If a window length is specified, split each query seq into windows of this size.
query_window_ids = []
query_window_seqs = []
if window_length != None:
    print("Splitting the query into windows...", end="")
    sys.stdout.flush()
    for x in range(len(query_ids)):
        window_start_pos = 0
        window_end_pos = window_length-1
        while window_end_pos < len(query_seqs[x]):
            window_pos_string = "_W" + str(window_start_pos) + "_" + str(window_end_pos)
            query_window_ids.append(query_ids[x] + window_pos_string)
            query_window_seqs.append(query_seqs[x][window_start_pos:window_end_pos+1])
            window_start_pos += window_shift
            window_end_pos += window_shift
    print(" done.")
else:
    query_window_ids = query_ids
    query_window_seqs = query_seqs

# Go through all sequence windows and use the sequence for BLAST searches against all subjects.
for z in range(len(query_window_ids)):

    # Determine the file name.
    if output_format == "phylip":
        aln_file_name = query_window_ids[z] + ".phy"
    else:
        aln_file_name = query_window_ids[z] + ".fasta"

    # If the file already exists, skip to next query window.
    if os.path.isfile(aln_file_name) and overwrite == False:
        print("Skipping blast search for file " + aln_file_name + " as file exists already.")
        continue

    # Prepare a string in fasta format for this query window.
    query_string = ">" + query_window_ids[z] + "\n" + query_window_seqs[z]

    # Write the sequence to a temporary file, so that it can be used as a query for BLAST searches.
    tmp_query = tempfile.NamedTemporaryFile(delete=False)
    tmp_query.write(query_string.encode('utf-8'))
    tmp_query.close()

    # Prepare lists to store the lists containing gap opening positions, and the lists
    # containing the aligned subject sequences.
    gaps_opened_in_query_per_subject = []
    sseqids_per_subject = []
    bitscores_per_subject = []
    if strictness > 0:
        number_of_hits_per_subject = []
    aligned_subject_seq_per_subject = []
    if translate:
        nucl_subject_seq_per_subject = []
        aligned_nucl_subject_seq_per_subject = []

    # Run BLAST searches against each subject.
    for subject in subjects:

        # Feedback.
        if window_length != None:
            print("Running blast searches for window " + query_window_ids[z] + " in subject " + subject + "...", end="")
            sys.stdout.flush()
        else:
            print("Running blast searches for query " + query_window_ids[z] + " in subject " + subject + "...", end="")
            sys.stdout.flush()

        # Prepare the BLAST command as array.
        FNULL = open(os.devnull, 'w')
        if translate:
            blast_ary = ["tblastn"]
        else:
            blast_ary = ["blastn"]
        if genetic_code != None:
            blast_ary.append("-db_gencode")
            blast_ary.append(str(genetic_code))
        blast_ary.append("-query")
        blast_ary.append(tmp_query.name)
        blast_ary.append("-db")
        blast_ary.append(subject)
        if query_evalues[z] != None:
            blast_ary.append("-evalue")
            blast_ary.append(str(query_evalues[z]))
        if strictness > 0:
            blast_ary.append("-culling_limit")
            blast_ary.append("1")
        blast_ary.append("-outfmt")
        blast_ary.append("6 qseqid sseqid evalue bitscore pident qstart qend sstart send qseq sseq")

        # Call BLAST.
        blast_out = subprocess.check_output(blast_ary).decode("utf-8").strip()
        blast_lines = blast_out.split("\n")

        if blast_lines == [""]:
            print("done.")
            aligned_subject_seq = ""
            for w in range(len(query_window_seqs[z])):
                aligned_subject_seq += "-"
            gaps_opened_in_query_per_subject.append([])
            if strictness == 0:
                sseqids_per_subject.append(["None"])
            else:
                sseqids_per_subject.append("None")
            if strictness == 0:
                bitscores_per_subject.append(["None"])
            else:
                bitscores_per_subject.append("None")
            if strictness > 0:
                number_of_hits_per_subject.append(0)
            aligned_subject_seq_per_subject.append(aligned_subject_seq)
            if translate:
                if strictness == 0:
                    nucl_subject_seq_per_subject.append([])
                else:
                    nucl_subject_seq_per_subject.append("")
        else:
            qseqids = []
            sseqids = []
            evalues = []
            bitscores = []
            pidents = []
            qstarts = []
            qends = []
            sstarts = []
            sends = []
            qseqs = []
            sseqs = []
            for blast_line in blast_lines:
                blast_line_ary = blast_line.split("\t")
                if len(blast_line_ary) == 11:
                    qseqids.append(blast_line_ary[0])
                    sseqids.append(blast_line_ary[1])
                    evalues.append(float(blast_line_ary[2]))
                    bitscores.append(float(blast_line_ary[3]))
                    pidents.append(float(blast_line_ary[4]))
                    qstarts.append(int(blast_line_ary[5]))
                    qends.append(int(blast_line_ary[6]))
                    sstarts.append(int(blast_line_ary[7]))
                    sends.append(int(blast_line_ary[8]))
                    qseqs.append(blast_line_ary[9])
                    sseqs.append(blast_line_ary[10])

            # If a bitscore threshold had been specified, remove all hits below
            # this threshold. (The same does not need to be done for evalue 
            # thresholds, as this can be filtered for when running the blast search).
            if query_bitscores[z] != None:
                tmp_qseqids = []
                tmp_sseqids = []
                tmp_evalues = []
                tmp_bitscores = []
                tmp_pidents = []
                tmp_qstarts = []
                tmp_qends = []
                tmp_sstarts = []
                tmp_sends = []
                tmp_qseqs = []
                tmp_sseqs = []
                for x in range(len(qseqids)):
                    if bitscores[x] >= query_bitscores[z]:
                        tmp_qseqids.append(qseqids[x])
                        tmp_sseqids.append(sseqids[x])
                        tmp_evalues.append(evalues[x])
                        tmp_bitscores.append(bitscores[x])
                        tmp_pidents.append(pidents[x])
                        tmp_qstarts.append(qstarts[x])
                        tmp_qends.append(qends[x])
                        tmp_sstarts.append(sstarts[x])
                        tmp_sends.append(sends[x])
                        tmp_qseqs.append(qseqs[x])
                        tmp_sseqs.append(sseqs[x])
                qseqids = tmp_qseqids
                sseqids = tmp_sseqids
                evalues = tmp_evalues
                bitscores = tmp_bitscores
                pidents = tmp_pidents
                qstarts = tmp_qstarts
                qends = tmp_qends
                sstarts = tmp_sstarts
                sends = tmp_sends
                qseqs = tmp_qseqs
                sseqs = tmp_sseqs

            # If tblastn was used, remove all hits that contain a stop codon (*).
            if translate:
                tmp_qseqids = []
                tmp_sseqids = []
                tmp_evalues = []
                tmp_bitscores = []
                tmp_pidents = []
                tmp_qstarts = []
                tmp_qends = []
                tmp_sstarts = []
                tmp_sends = []
                tmp_qseqs = []
                tmp_sseqs = []
                for x in range(len(qseqids)):
                    if "*" not in sseqs[x]:
                        tmp_qseqids.append(qseqids[x])
                        tmp_sseqids.append(sseqids[x])
                        tmp_evalues.append(evalues[x])
                        tmp_bitscores.append(bitscores[x])
                        tmp_pidents.append(pidents[x])
                        tmp_qstarts.append(qstarts[x])
                        tmp_qends.append(qends[x])
                        tmp_sstarts.append(sstarts[x])
                        tmp_sends.append(sends[x])
                        tmp_qseqs.append(qseqs[x])
                        tmp_sseqs.append(sseqs[x])
                qseqids = tmp_qseqids
                sseqids = tmp_sseqids
                evalues = tmp_evalues
                bitscores = tmp_bitscores
                pidents = tmp_pidents
                qstarts = tmp_qstarts
                qends = tmp_qends
                sstarts = tmp_sstarts
                sends = tmp_sends
                qseqs = tmp_qseqs
                sseqs = tmp_sseqs

            # If not a single blast line was added (after bitscore filtering),
            # return an empty sequence.
            if qseqids == []:
                aligned_subject_seq = ""
                for w in range(len(query_window_seqs[z])):
                    aligned_subject_seq += "-"
                gaps_opened_in_query_per_subject.append([])
                if strictness == 0:
                    sseqids_per_subject.append(["None"])
                else:
                    sseqids_per_subject.append("None")
                if strictness == 0:
                    bitscores_per_subject.append(["None"])
                else:
                    bitscores_per_subject.append("None")
                if strictness > 0:
                    number_of_hits_per_subject.append(0)
                aligned_subject_seq_per_subject.append(aligned_subject_seq)
                if translate:
                    if strictness == 0:
                        nucl_subject_seq_per_subject.append([])
                    else:
                        nucl_subject_seq_per_subject.append("")
                print(" done.")

            # Otherwise continue analysing the blast hits.
            else:

                # Find the contig/scaffold with most hits (only with strictness = 2).
                if strictness == 2:
                    uniq_sseqids = sorted(list(set(sseqids)))
                    uniq_sseqid_occurrences = []
                    for uniq_sseqid in uniq_sseqids:
                        uniq_sseqid_occurrences.append(sseqids.count(uniq_sseqid))
                    best_sseqid = uniq_sseqids[uniq_sseqid_occurrences.index(max(uniq_sseqid_occurrences))]

                # Remove hits for all other contigs/scaffolds (only with strictness = 2).
                if strictness == 2:
                    tmp_qseqids = []
                    tmp_sseqids = []
                    tmp_evalues = []
                    tmp_bitscores = []
                    tmp_pidents = []
                    tmp_qstarts = []
                    tmp_qends = []
                    tmp_sstarts = []
                    tmp_sends = []
                    tmp_qseqs = []
                    tmp_sseqs = []
                    for x in range(len(qseqids)):
                        if sseqids[x] == best_sseqid:
                            tmp_qseqids.append(qseqids[x])
                            tmp_sseqids.append(sseqids[x])
                            tmp_evalues.append(evalues[x])
                            tmp_bitscores.append(bitscores[x])
                            tmp_pidents.append(pidents[x])
                            tmp_qstarts.append(qstarts[x])
                            tmp_qends.append(qends[x])
                            tmp_sstarts.append(sstarts[x])
                            tmp_sends.append(sends[x])
                            tmp_qseqs.append(qseqs[x])
                            tmp_sseqs.append(sseqs[x])
                    qseqids = tmp_qseqids
                    sseqids = tmp_sseqids
                    evalues = tmp_evalues
                    bitscores = tmp_bitscores
                    pidents = tmp_pidents
                    qstarts = tmp_qstarts
                    qends = tmp_qends
                    sstarts = tmp_sstarts
                    sends = tmp_sends
                    qseqs = tmp_qseqs
                    sseqs = tmp_sseqs

                # Determine the orientation for all hits (only with strictness = 2).
                if strictness == 2:
                    sorients = []
                    for x in range(len(qseqids)):
                        if qstarts[x] < qends[x]:
                            if sstarts[x] < sends[x]:
                                sorients.append("for")
                            else:
                                sorients.append("rev")
                        else:
                            if sstarts[x] < sends[x]:
                                sorients.append("rev")
                            else:
                                sorients.append("for")

                # Find the more common orientation (only with strictness = 2).
                if strictness == 2:
                    if sorients.count("rev") > sorients.count("for"):
                        common_orient = "rev"
                    else:
                        common_orient = "for"

                # Discard hits that have the less common orientation (only with strictness = 2).
                if strictness == 2:
                    tmp_qseqids = []
                    tmp_sseqids = []
                    tmp_evalues = []
                    tmp_bitscores = []
                    tmp_pidents = []
                    tmp_qstarts = []
                    tmp_qends = []
                    tmp_sstarts = []
                    tmp_sends = []
                    tmp_qseqs = []
                    tmp_sseqs = []
                    for x in range(len(qseqids)):
                        if sorients[x] == common_orient:
                            tmp_qseqids.append(qseqids[x])
                            tmp_sseqids.append(sseqids[x])
                            tmp_evalues.append(evalues[x])
                            tmp_bitscores.append(bitscores[x])
                            tmp_pidents.append(pidents[x])
                            tmp_qstarts.append(qstarts[x])
                            tmp_qends.append(qends[x])
                            tmp_sstarts.append(sstarts[x])
                            tmp_sends.append(sends[x])
                            tmp_qseqs.append(qseqs[x])
                            tmp_sseqs.append(sseqs[x])
                    qseqids = tmp_qseqids
                    sseqids = tmp_sseqids
                    evalues = tmp_evalues
                    bitscores = tmp_bitscores
                    pidents = tmp_pidents
                    qstarts = tmp_qstarts
                    qends = tmp_qends
                    sstarts = tmp_sstarts
                    sends = tmp_sends
                    qseqs = tmp_qseqs
                    sseqs = tmp_sseqs

                # Sort all arrays by the qstart.
                all_sorted = False
                while all_sorted == False:
                    all_sorted = True
                    for x in range(len(qseqids)-1):
                        if qstarts[x] > qstarts[x+1]:
                            all_sorted = False
                            qseqids[x], qseqids[x+1] = qseqids[x+1], qseqids[x]
                            sseqids[x], sseqids[x+1] = sseqids[x+1], sseqids[x]
                            evalues[x], evalues[x+1] = evalues[x+1], evalues[x]
                            bitscores[x], bitscores[x+1] = bitscores[x+1], bitscores[x]
                            pidents[x], pidents[x+1] = pidents[x+1], pidents[x]
                            qstarts[x], qstarts[x+1] = qstarts[x+1], qstarts[x]
                            qends[x], qends[x+1] = qends[x+1], qends[x]
                            sstarts[x], sstarts[x+1] = sstarts[x+1], sstarts[x]
                            sends[x], sends[x+1] = sends[x+1], sends[x]
                            qseqs[x], qseqs[x+1] = qseqs[x+1], qseqs[x]
                            sseqs[x], sseqs[x+1] = sseqs[x+1], sseqs[x]

                # If two or more hits have overlapping query ranges, discard the one with lower evalue.
                if strictness > 0:
                    no_overlaps = False
                    while no_overlaps == False:
                        no_overlaps = True
                        for x in range(len(qseqids)-1):
                            if qends[x] >= qstarts[x+1]:
                                no_overlaps = False
                                if evalues[x] < evalues[x+1]:
                                    # Remove hit x+1 from the dataset.
                                    del qseqids[x+1]
                                    del sseqids[x+1]
                                    del evalues[x+1]
                                    del bitscores[x+1]
                                    del pidents[x+1]
                                    del qstarts[x+1]
                                    del qends[x+1]
                                    del sstarts[x+1]
                                    del sends[x+1]
                                    del qseqs[x+1]
                                    del sseqs[x+1]
                                else:
                                    # Remove hit x from the dataset.
                                    del qseqids[x]
                                    del sseqids[x]
                                    del evalues[x]
                                    del bitscores[x]
                                    del pidents[x]
                                    del qstarts[x]
                                    del qends[x]
                                    del sstarts[x]
                                    del sends[x]
                                    del qseqs[x]
                                    del sseqs[x]
                                break

                # Determine the order for all queries (only with strictness = 2).
                if strictness == 2:
                    qorder = []
                    order = 1
                    for x in range(len(qseqids)):
                        qorder.append(order)
                        order += 1

                # Determine the order for all subjects (only with strictness = 2).
                if strictness == 2:
                    sorder = []
                    sorted_sstarts = sorted(sstarts)
                    if common_orient == "rev":
                        sorted_sstarts.reverse()
                    for x in range(len(qseqids)):
                        sorder.append(sorted_sstarts.index(sstarts[x])+1)

                # Discard those hits that violate the order until the two orders are the same
                # (only with strictness = 2).
                if strictness == 2:
                    while sorder != qorder:

                        # Find the hits that violate the order the most.
                        violation = []
                        for x in range(len(qseqids)):
                            violation.append(abs(sorder[x]-qorder[x]))
                        max_violation = max(violation)

                        # Discard hits that violate the order the most.
                        tmp_qseqids = []
                        tmp_sseqids = []
                        tmp_evalues = []
                        tmp_bitscores = []
                        tmp_pidents = []
                        tmp_qstarts = []
                        tmp_qends = []
                        tmp_sstarts = []
                        tmp_sends = []
                        tmp_qseqs = []
                        tmp_sseqs = []
                        for x in range(len(qseqids)):
                            if violation[x] < max_violation:
                                tmp_qseqids.append(qseqids[x])
                                tmp_sseqids.append(sseqids[x])
                                tmp_evalues.append(evalues[x])
                                tmp_bitscores.append(bitscores[x])
                                tmp_pidents.append(pidents[x])
                                tmp_qstarts.append(qstarts[x])
                                tmp_qends.append(qends[x])
                                tmp_sstarts.append(sstarts[x])
                                tmp_sends.append(sends[x])
                                tmp_qseqs.append(qseqs[x])
                                tmp_sseqs.append(sseqs[x])
                        qseqids = tmp_qseqids
                        sseqids = tmp_sseqids
                        evalues = tmp_evalues
                        bitscores = tmp_bitscores
                        pidents = tmp_pidents
                        qstarts = tmp_qstarts
                        qends = tmp_qends
                        sstarts = tmp_sstarts
                        sends = tmp_sends
                        qseqs = tmp_qseqs
                        sseqs = tmp_sseqs

                        # Determine the order for all queries.
                        qorder = []
                        order = 1
                        for x in range(len(qseqids)):
                            qorder.append(order)
                            order += 1

                        # Determine the order for all subjects.
                        sorder = []
                        sorted_sstarts = sorted(sstarts)
                        if common_orient == "rev":
                            sorted_sstarts.reverse()
                        for x in range(len(qseqids)):
                            sorder.append(sorted_sstarts.index(sstarts[x])+1)

                # Initiate gaps_opened_in_query, a list to store the position of all gaps that
                # this subject (with strictness > 0) / these subjects (with strictness = 0)
                # cause(s) in the query sequence.
                # With strictness = 0, this list contains other lists, one for each
                # subject. 
                if strictness == 0:

                    # Initiate gaps_opened_in_query
                    gaps_opened_in_query = []

                    # Initiate a list to store the starting positions of each subject sequence in
                    # a multiple sequence alignment with gaps.
                    sstarts_in_msa = []
                    
                    # In order to know which gaps are opened in the query by each subject sequence,
                    # the query sequences need to be checked for "-".
                    for x in range(len(qseqids)):
                        
                        # Make the xth item of list gaps_opened_in_query a list in itself.
                        gaps_opened_in_query.append([])
                        
                        # Find out how many gaps have already been opened before the start position
                        # of this query sequence region (=subject sequence region).
                        gaps_opened_before_query_start = []
                        for alist in gaps_opened_in_query:
                            for item in alist:
                                if item < qstarts[x]:
                                    if item not in gaps_opened_before_query_start:
                                        gaps_opened_before_query_start.append(item)
                        number_of_gaps_opened_before_query_start = len(gaps_opened_before_query_start)

                        # For this query/subject sequence region, determine its start position in a
                        # multiple sequence alignment containing gaps.
                        sstarts_in_msa.append(qstarts[x] + number_of_gaps_opened_before_query_start)

                        # For each position in the query/subject sequence region, see whether the subject
                        # has opened a gap in its query. If so, if it's already stored in gaps_opened_in_query,
                        # only register that this subject is also responsible for it (by storing it in
                        # the list of gaps_opened_in_query that corresponds to this subject). If it's not yet
                        # stored in gaps_opened_in_query, do so and shift all gaps at positions greater than
                        # this one up by one.
                        number_of_gaps_opened_up_to_here = number_of_gaps_opened_before_query_start
                        for y in range(len(qseqs[x])):
                            pos_in_region = y
                            pos_in_alignment = y + qstarts[x] + number_of_gaps_opened_up_to_here
                            # See whether there is a gap already known at this position in the alignment.
                            gap_known_at_this_position = False
                            for alist in gaps_opened_in_query:
                                if pos_in_alignment in alist:
                                    gap_known_at_this_position = True
                            # See whether there is a gap found at this position in this query:
                            gap_found_in_this_query_at_this_position = False
                            if qseqs[x][y] == "-":
                                gap_found_in_this_query_at_this_position = True
                            # 1.) If there is a gap already registered at this position and it is also found in this
                            #     query, then memorize that this subject is also responsible for the observed
                            #     gap (don't change number_of_gaps_opened_up_to_here).
                            if gap_known_at_this_position == True and gap_found_in_this_query_at_this_position == True:
                                gaps_opened_in_query[x].append(pos_in_alignment)
                            # 2.) If there is no gap known for this position, but there is one found in this query,
                            #     then memorize that this subject is responsible for the gap, and shift all gaps
                            #     at positions greater than this one up by one (don't change number_of_gaps_opened_up_to_here).
                            elif gap_known_at_this_position == False and gap_found_in_this_query_at_this_position == True:
                                gaps_opened_in_query[x].append(pos_in_alignment)
                                for alist in gaps_opened_in_query:
                                    for item in alist:
                                        if item > pos_in_alignment:
                                            item += 1
                            # 3.) If there is a gap known for this position, but this query does not have one,
                            #     increase the number_of_gaps_opened_up_to_here by one.
                            elif gap_known_at_this_position == True and gap_found_in_this_query_at_this_position == False:
                                number_of_gaps_opened_up_to_here += 1
                            # 4.) If there is no gap known for this position, and this query also does not have one,
                            #     move on.
                else:
                    # If all subjects are sequential, as they are with strictness > 0, it's much easier.
                    gaps_opened_in_query = []
                    sstarts_in_msa = []
                    for x in range(len(qseqids)):
                        number_of_gaps_already_opened = len(gaps_opened_in_query)
                        sstarts_in_msa.append(qstarts[x] + number_of_gaps_already_opened)
                        for y in range(len(qseqs[x])):
                            if qseqs[x][y] == "-":
                                gaps_opened_in_query.append(qstarts[x]+y+number_of_gaps_already_opened)

                # Prepare the aligned sequence for this subject.
                if strictness == 0:
                    aligned_subject_seqs = []
                    for x in range(len(qseqids)):
                        aligned_subject_seq = ""
                        while len(aligned_subject_seq) < sstarts_in_msa[x]-1:
                            aligned_subject_seq += "-"
                        pos_in_region = 0
                        pos_in_alignment = sstarts_in_msa[x]
                        while pos_in_region < len(sseqs[x]):
                            # See whether any subject has opened a gap in the alignment at this position.
                            gap_known_at_this_position = False
                            for alist in gaps_opened_in_query:
                                if pos_in_alignment in alist:
                                    gap_known_at_this_position = True
                            # See whether this subject is responsible for the gap in the alignment at this position.
                            this_subject_is_responsible = False
                            if pos_in_alignment in gaps_opened_in_query[x]:
                                this_subject_is_responsible = True
                            # 1.) If there is a gap in the alignment at this position and this subject
                            #     is responsible for it, just add this position of this subject to the
                            #     aligned sequence and move on.
                            if gap_known_at_this_position == True and this_subject_is_responsible == True:
                                aligned_subject_seq += sseqs[x][pos_in_region]
                                pos_in_region += 1
                                pos_in_alignment += 1
                            # 2.) If there is a gap in the alignment at this position and this subject is
                            #     not responsible for it, add a gap to the aligned sequence.
                            elif gap_known_at_this_position == True and this_subject_is_responsible == False:
                                aligned_subject_seq += "-"
                                pos_in_alignment += 1
                            # 3.) If there is no gap in the alignment at this position, just add this position
                            #     of this subject to the aligned sequence and move on.
                            elif gap_known_at_this_position == False:
                                aligned_subject_seq += sseqs[x][pos_in_region]
                                pos_in_region += 1
                                pos_in_alignment += 1
                        
                        # See how many gaps have been opened in total.
                        gaps_opened_in_query_flat = []
                        for alist in gaps_opened_in_query:
                            for item in alist:
                                if item not in gaps_opened_in_query_flat:
                                    gaps_opened_in_query_flat.append(item)
                        total_number_of_opened_gaps = len(gaps_opened_in_query_flat)
                        alignment_length = len(query_window_seqs[z]) + total_number_of_opened_gaps

                        # Add missing data to the end of the alignment.
                        while pos_in_alignment <= alignment_length:
                            aligned_subject_seq += "-"
                            pos_in_alignment += 1

                        # Store the aligned sequence for this subject's sequence.
                        aligned_subject_seqs.append(aligned_subject_seq)

                else:
                    aligned_subject_seq = ""
                    for x in range(len(qseqids)):
                        while len(aligned_subject_seq) < sstarts_in_msa[x]-1:
                            aligned_subject_seq += "-"
                        aligned_subject_seq += sseqs[x]
                    while len(aligned_subject_seq) < len(query_window_seqs[z]) + len(gaps_opened_in_query):
                        aligned_subject_seq += "-"

                # Combine the subject ids in a string to be included in the fasta id.
                if strictness > 0:
                    if strictness == 1:
                        combined_sseqids = ""
                        for sseqid in sseqids:
                            combined_sseqids += sseqid
                            combined_sseqids += "|"
                        combined_sseqids = combined_sseqids[:-1]
                    elif strictness == 2:
                        combined_sseqids == sseqid[0]
                    else:
                        print("ERROR: Unknown strictness setting: " + strictness + "!")
                        sys.exit(1)

                # Calculate the combined bitscore to be included in the fasta id.
                if strictness > 0:
                    combined_bitscore = 0
                    for bitscore in bitscores:
                        combined_bitscore += bitscore

                # Store the results from this subject.
                # sseqids_per_subject, bitscores_per_subject, aligned_subject_seq_per_subject all represent
                # lists of lists when strictness is 0, and lists of values when strictness is 1 or 2.
                if strictness == 0:
                    gaps_opened_in_query_flat = []
                    for alist in gaps_opened_in_query:
                        for item in alist:
                            if item not in gaps_opened_in_query_flat:
                                gaps_opened_in_query_flat.append(item)
                    gaps_opened_in_query_per_subject.append(gaps_opened_in_query_flat)
                else:
                    gaps_opened_in_query_per_subject.append(gaps_opened_in_query)
                if strictness == 0:
                    sseqids_per_subject.append(sseqids)
                else:
                    sseqids_per_subject.append(combined_sseqids)
                if strictness == 0:
                    bitscores_per_subject.append(bitscores)
                else:
                    bitscores_per_subject.append(combined_bitscore)
                if strictness > 0:
                    number_of_hits_per_subject.append(len(qseqids))
                if strictness == 0:
                    aligned_subject_seq_per_subject.append(aligned_subject_seqs)
                else:
                    aligned_subject_seq_per_subject.append(aligned_subject_seq)

                # If tblastn was used, the nucleotide sequence has to be obtained from the subject
                # fasta file.
                if translate:
                    print(" done.")
                    print("Reading subject " + subject + " to retrieve nucleotide sequence...", end="")
                    sys.stdout.flush()
                    if strictness == 0:
                        nucl_subject_seqs = []
                    else:
                        nucl_subject_seq = ""
                    for x in range(len(qseqids)):
                        with open(subject) as f:
                            subject_nucl_seq = ""
                            in_seq = False
                            for line in f:
                                if line[0:len(sseqids[x])+1] == ">" + sseqids[x]:
                                    if line[len(sseqids[x])+1] == " " or line[len(sseqids[x])+1] == "\n":
                                        in_seq = True
                                elif in_seq == True and line.strip() != "":
                                    if line[0] == ">":
                                        break
                                    else:
                                        subject_nucl_seq += line.strip()
                            if subject_nucl_seq == "":
                                print("ERROR: The nucleotide sequence " + sseqids[x] + " could not be found in subject " + subject + "!")
                                sys.exit(1)
                            if sstarts[x] < sends[x]:
                                nucl_subject_seq_for_this_hit = subject_nucl_seq[sstarts[x]-1:sends[x]]
                            else:
                                nucl_subject_seq_for_this_hit = subject_nucl_seq[sends[x]-1:sstarts[x]]
                                # Make reverse complement.
                                tmp = ""
                                pos = len(nucl_subject_seq_for_this_hit)
                                while pos > 0:
                                    pos -= 1
                                    if nucl_subject_seq_for_this_hit[pos].upper() == "A":
                                        tmp += "T"
                                    elif nucl_subject_seq_for_this_hit[pos].upper() == "C":
                                        tmp += "G"
                                    elif nucl_subject_seq_for_this_hit[pos].upper() == "G":
                                        tmp += "C"
                                    elif nucl_subject_seq_for_this_hit[pos].upper() == "T":
                                        tmp += "A"
                                    elif nucl_subject_seq_for_this_hit[pos].strip() != "":
                                        tmp += "N"
                                nucl_subject_seq_for_this_hit = tmp
                            if strictness == 0:
                                nucl_subject_seqs.append(nucl_subject_seq_for_this_hit)
                            else:
                                nucl_subject_seq += nucl_subject_seq_for_this_hit
                    if strictness == 0:
                        nucl_subject_seq_per_subject.append(nucl_subject_seqs)
                    else:
                        nucl_subject_seq_per_subject.append(nucl_subject_seq)
                print(" done.")

    # Produce the alignment with the query and all subject sequences.
    print("Producing alignment...", end="")
    reached_the_end = False
    pos = 0
    aligned_ref_seq = query_window_seqs[z]
    while reached_the_end == False:
        # Increment the current position.
        pos += 1
        # Test whether any (and if so, which) of the subject has opened a gap at this position.
        subjects_responsible_for_gap_opening_at_this_pos = []
        for subject_count in range(len(subjects)):
            if pos in gaps_opened_in_query_per_subject[subject_count]:
                subjects_responsible_for_gap_opening_at_this_pos.append(subject_count)

        # If there is indeed a gap opening at this position.
        if len(subjects_responsible_for_gap_opening_at_this_pos) > 0:

            # Insert a gap into the query sequence.
            tmp_aligned_ref_seq = aligned_ref_seq[:(pos-1)] + "-" + aligned_ref_seq[(pos-1):]
            aligned_ref_seq = tmp_aligned_ref_seq

            # Insert a gap into all subject sequences unless they are responsible for the gap.
            for subject_count in range(len(subjects)):
                if subject_count not in subjects_responsible_for_gap_opening_at_this_pos:
                    if strictness == 0:
                        for sseq_count in range(len(aligned_subject_seq_per_subject[subject_count])):
                            aligned_subject_seq in aligned_subject_seq_per_subject[subject_count][sseq_count]
                            tmp_aligned_subject_seq = aligned_subject_seq[:(pos-1)] + "-" + aligned_subject_seq[(pos-1):]
                            aligned_subject_seq_per_subject[subject_count][sseq_count] = tmp_aligned_subject_seq
                            # Increase all gap positions greater than the current pos by 1.
                            tmp_gaps_opened_in_query = []
                            for gap_pos in gaps_opened_in_query_per_subject[subject_count]:
                                if gap_pos > pos:
                                    tmp_gaps_opened_in_query.append(gap_pos+1)
                                else:
                                    tmp_gaps_opened_in_query.append(gap_pos)
                            gaps_opened_in_query_per_subject[subject_count] = tmp_gaps_opened_in_query
                    else:
                        aligned_subject_seq = aligned_subject_seq_per_subject[subject_count]
                        tmp_aligned_subject_seq = aligned_subject_seq[:(pos-1)] + "-" + aligned_subject_seq[(pos-1):]
                        aligned_subject_seq_per_subject[subject_count] = tmp_aligned_subject_seq
                        # Increase all gap positions greater than the current pos by 1.
                        tmp_gaps_opened_in_query = []
                        for gap_pos in gaps_opened_in_query_per_subject[subject_count]:
                            if gap_pos > pos:
                                tmp_gaps_opened_in_query.append(gap_pos+1)
                            else:
                                tmp_gaps_opened_in_query.append(gap_pos)
                        gaps_opened_in_query_per_subject[subject_count] = tmp_gaps_opened_in_query

        # Stop when the current position has reached the length of the aligned query seq.
        if pos >= len(aligned_ref_seq):
            reached_the_end = True

    # Collect the alignment entries.
    aln_ids = []
    aln_seqs = []
    aln_ids.append(query_window_ids[z])
    aln_seqs.append(aligned_ref_seq)
    if strictness == 0:
        for subject_count in range(len(subjects)):
            for sseqids_count in range(len(sseqids_per_subject[subject_count])):
                aln_id = subjects[subject_count].split("/")[-1].replace(".fasta","").replace(".ctg","").replace(".utg","")
                aln_id += "[&sseqid=" + str(sseqids_per_subject[subject_count][sseqids_count]) + ","
                aln_id += "bitscore=" + str(bitscores_per_subject[subject_count][sseqids_count])
                aln_id += "]"
                aln_ids.append(aln_id)
                aln_seqs.append(aligned_subject_seq_per_subject[subject_count][sseqids_count])
    else:
        for subject_count in range(len(subjects)):
            aln_id = subjects[subject_count].split("/")[-1].replace(".fasta","").replace(".ctg","").replace(".utg","")
            aln_id += "[&sseqid=" + str(sseqids_per_subject[subject_count]) + ","
            aln_id += "bitscore=" + str(bitscores_per_subject[subject_count]) + ","
            aln_id += "nhits=" + str(number_of_hits_per_subject[subject_count])
            aln_id += "]"
            aln_ids.append(aln_id)
            aln_seqs.append(aligned_subject_seq_per_subject[subject_count])
    print(" done.")

    if refine:
        # User mafft to refine the alignment.
        print("Refining the alignment with mafft...", end="")
        tmp_aln_string = ""
        for aln_count in range(len(aln_ids)):
            tmp_aln_string += ">" + aln_ids[aln_count].replace(" ","_")
            tmp_aln_string += "\n"
            tmp_aln_string += aln_seqs[aln_count]
            tmp_aln_string += "\n"
        tmp_alignment_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_alignment_file.write(tmp_aln_string.encode('utf-8'))
        tmp_alignment_file.close()
        mafft_ary = ["mafft"]
        mafft_ary.append("--auto")
        mafft_ary.append("--quiet")
        mafft_ary.append(tmp_alignment_file.name)
        mafft_out = subprocess.check_output(mafft_ary).decode("utf-8").strip()
        refined_aln_ids = []
        refined_aln_seqs = []
        mafft_lines = mafft_out.split("\n")
        for mafft_line in mafft_lines:
            if mafft_line[0] == ">":
                refined_aln_ids.append(mafft_line.strip()[1:])
                refined_aln_seqs.append("")
            elif mafft_line.strip() != "":
                refined_aln_seqs[-1] += mafft_line.strip()
        # Make sure that the alignment ids haven't changed.
        if aln_ids != refined_aln_ids:
            print("ERROR: Alignment ids have changed after mafft refinement! (before: " + aln_ids[0] + ", " + aln_ids[1] + ", " + aln_ids[2] + "..." + "; after:" + refined_aln_ids[0] + ", " + refined_aln_ids[1] + ", " + refined_aln_ids[2] + "...)")
            sys.exit(1)
        if len(aln_seqs) != len(refined_aln_seqs):
            print("ERROR: The number of sequences has changed after mafft refinement!")
        for x in range(len(aln_seqs)):
            if aln_seqs[x].replace("-","") != refined_aln_seqs[x].replace("-",""):
                print("ERROR: The composition of sequence " + aln_ids[x] + " has changed after mafft refinement!")
                print("Before mafft refinement:")
                print(aln_seqs[x])
                print()
                print("After mafft refinement:")
                print(refined_aln_seqs[x])
                sys.exit(1)
        aln_ids = refined_aln_ids
        aln_seqs = refined_aln_seqs
        print(" done.")

    # Check whether the minimum completeness is met by this alignment.
    write_this_alignment = True
    if minimum_completeness != None:
        print("Checking whether alignment completeness is sufficient...", end="")
        number_of_complete_sites = 0
        for pos in range(len(aln_seqs[0])):
            # Check whether there are gaps at this position.
            gaps_at_this_pos = False
            for aln_seq in aln_seqs:
                if aln_seq[pos] == "-":
                    gaps_at_this_pos = True
                    break
            # If there are no gaps at this position, increase the counter.
            if gaps_at_this_pos == False:
                number_of_complete_sites += 1
        if number_of_complete_sites < minimum_completeness:
            write_this_alignment = False
        print(" done.")

    # If tblastn was used, use the amino acid alignment and the obtained nucleotide sequence
    # to produce aligned nucleotide sequences. Leave out the first sequence, which is the
    # amino acid query.
    if translate and write_this_alignment:
        print("Producing nucleotide alignment based on amino acid alignment...", end="")
        untranslated_aln_ids = aln_ids[1:]
        untranslated_aln_seqs = []
        aln_index = 0
        for x in range(len(nucl_subject_seq_per_subject)):
            if strictness == 0:
                for y in range(len(nucl_subject_seq_per_subject[x])):
                    aln_index += 1
                    aligned_nucl_subject_seq_per_subject = ""
                    # Make sure that the nucleotide sequence contains three times the number of sites
                    # of the amino acid sequence.
                    if len(nucl_subject_seq_per_subject[x][y]) != 3 * len(aln_seqs[aln_index].replace("-","")):
                        print("ERROR: The length of the amino acid and nucleotide sequcences don't match!")
                        print("ERROR: " + aln_ids[aln_index])
                        print("ERROR: " + nucl_subject_seq_per_subject[x][y])
                        print("ERROR: " + aln_seqs[aln_index])
                        sys.exit(1)
                    nucl_pos = 0
                    for aa_pos in range(len(aln_seqs[aln_index])):
                        if aln_seqs[aln_index][aa_pos] == "-":
                            aligned_nucl_subject_seq_per_subject += "---"
                        else:
                            aligned_nucl_subject_seq_per_subject += nucl_subject_seq_per_subject[x][y][nucl_pos:nucl_pos+3]
                            nucl_pos += 3
                    untranslated_aln_seqs.append(aligned_nucl_subject_seq_per_subject)
            else:
                aln_index += 1
                aligned_nucl_subject_seq_per_subject = ""
                # Make sure that the nucleotide sequence contains three times the number of sites
                # of the amino acid sequence.
                if len(nucl_subject_seq_per_subject[x]) != 3 * len(aln_seqs[aln_index].replace("-","")):
                    print("ERROR: The length of the amino acid and nucleotide sequcences don't match!")
                    print("ERROR: " + aln_ids[aln_index])
                    print("ERROR: " + nucl_subject_seq_per_subject[x])
                    print("ERROR: " + aln_seqs[aln_index])
                    sys.exit(1)
                nucl_pos = 0
                for aa_pos in range(len(aln_seqs[x+1])):
                    if aln_seqs[x+1][aa_pos] == "-":
                        aligned_nucl_subject_seq_per_subject += "---"
                    else:
                        aligned_nucl_subject_seq_per_subject += nucl_subject_seq_per_subject[x][nucl_pos:nucl_pos+3]
                        nucl_pos += 3
                untranslated_aln_seqs.append(aligned_nucl_subject_seq_per_subject)
        # Remove gap-only sites.
        if len(untranslated_aln_seqs) > 0:
            tmp_untranslated_aln_seqs = []
            for untranslated_aln_seq in untranslated_aln_seqs:
                tmp_untranslated_aln_seqs.append("")
            for x in range(len(untranslated_aln_seqs[0])):
                gap_only_site = True
                for untranslated_aln_seq in untranslated_aln_seqs:
                    if untranslated_aln_seq[x] != "-":
                        gap_only_site = False
                        break
                if gap_only_site == False:
                    for y in range(len(untranslated_aln_seqs)):
                        tmp_untranslated_aln_seqs[y] += untranslated_aln_seqs[y][x]
            untranslated_aln_seqs = tmp_untranslated_aln_seqs
        print(" done.")

    # If the alignment is ok, prepare it.
    aln_string = ""
    if translate:
        untranslated_aln_string = ""
    if write_this_alignment == True:

        # Prepare the alignment string.
        if output_format == "phylip":
            aln_string += str(len(aln_ids)) + " " + str(len(aln_seqs[0])) + "\n"
            longest_id_length = 0
            for aln_id in aln_ids:
                if len(aln_id) > longest_id_length:
                    longest_id_length = len(aln_id)
            for aln_count in range(len(aln_ids)):
                aln_string += aln_ids[aln_count].replace(" ","_").ljust(longest_id_length+2)
                aln_string += aln_seqs[aln_count]
                aln_string += "\n"
            # If tblastn was used, also prepare the nucleotide sequence alignment, which has one
            # record less than the amino acid sequene alignment, as the amino acid query is not
            # included.
            if translate:
                if len(untranslated_aln_seqs) > 0:
                    untranslated_aln_string += str(len(untranslated_aln_ids)) + " " + str(len(untranslated_aln_seqs[0])) + "\n"
                    for untranslated_aln_count in range(len(untranslated_aln_ids)):
                        untranslated_aln_string += untranslated_aln_ids[untranslated_aln_count].replace(" ","_").ljust(longest_id_length+2)
                        untranslated_aln_string += untranslated_aln_seqs[untranslated_aln_count]
                        untranslated_aln_string += "\n"

        else:
            for aln_count in range(len(aln_ids)):
                aln_string += ">" + aln_ids[aln_count].replace(" ","_")
                aln_string += "\n"
                aln_string += aln_seqs[aln_count]
                aln_string += "\n"
            # If tblastn was used, also prepare the nucleotide sequence alignment, which has one
            # record less than the amino acid sequene alignment, as the amino acid query is not
            # included.
            if translate:
                if len(untranslated_aln_seqs) > 0:
                    for untranslated_aln_count in range(len(untranslated_aln_ids)):
                        untranslated_aln_string += ">" + untranslated_aln_ids[untranslated_aln_count].replace(" ","_")
                        untranslated_aln_string += "\n"
                        untranslated_aln_string += untranslated_aln_seqs[untranslated_aln_count]
                        untranslated_aln_string += "\n"

    # Write the alignment to file.
    if write_this_alignment or write_empty_files:
        aln_file = open(aln_file_name, 'w')
        aln_file.write(aln_string)
        aln_file.close()
        # If tblastn was used, also write the nucleotide sequence alignment file.
        if translate:
            untranslated_aln_file_name = aln_file_name.replace(".fasta","_nucl.fasta")
            untranslated_aln_file = open(untranslated_aln_file_name, 'w')
            untranslated_aln_file.write(untranslated_aln_string)
            untranslated_aln_file.close()

    # Feedback.
    if write_this_alignment:
        print("Wrote file " + aln_file_name)
        if translate:
            print("Wrote file " + untranslated_aln_file_name)
        sys.stdout.flush()            
    else:
        print("Insufficient completeness (" + str(number_of_complete_sites) + " sites) for " + query_window_ids[z])
        sys.stdout.flush()
        if write_empty_files:
            print("Wrote empty file " + aln_file_name)
            if translate:
                print("Wrote empty file " + untranslated_aln_file_name)
            sys.stdout.flush()            
