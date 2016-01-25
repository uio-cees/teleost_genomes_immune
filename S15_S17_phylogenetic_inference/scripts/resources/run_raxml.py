#!/usr/local/bin/python3

# Michael Matschiner, 2015-07-09
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse, textwrap, random, tempfile, os
from subprocess import call
from Bio import AlignIO

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      Wrapper for raxml that allows piping.
      Produces the ML tree as a simple string, or included in
      a nexus format string with the sequence alignment. 
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.996'
    )
parser.add_argument(
    '-n', '--nexus',
    action='store_true',
    help='produce a nexus string including the ML tree'
    )
parser.add_argument(
    '-p', '--parsimony',
    action='store_true',
    dest='parsimony',
    help='use the RAxML parsimony tree instead of the ML tree (RAxML option -y)'
    )
parser.add_argument(
    '--one-category',
    action='store_true',
    dest='one_category',
    help='use a single rate category (RAxML option -m GTRCAT -V)'
    )
parser.add_argument(
    '-f', '--from',
    nargs=1,
    type=int,
    default=[1],
    dest='start',
    help="start position of analysis window (first = 1)"
    )
parser.add_argument(
    '-t', '--to',
    nargs=1,
    type=int,
    default=[-1],
    dest='end',
    help="end position of analysis window (first = 1)"
    )
parser.add_argument(
    '-x', '--exclude-taxon',
    nargs=1,
    type=str,
    default=[None],
    dest='exclude_taxon',
    help="exclude all taxa that match this name"
    )
parser.add_argument(
    '-o', '--outgroup',
    nargs=1,
    type=str,
    default=[None],
    dest='outgroup',
    help="use this taxon as the outgroup"
    )
parser.add_argument(
    '-q', '--partitions',
    nargs=1,
    type=str,
    default=[None],
    dest='partitions',
    help="specify the name of a partitions file"
    )
parser.add_argument(
    '-c', '--constraint',
    nargs=1,
    type=str,
    default=[None],
    dest='constraint',
    help="specify the name of a constraint file"
    )
parser.add_argument(
    '-b', '--bootstraps',
    nargs=1,
    type=str,
    default=[None],
    dest='bootstraps',
    help="assess node support with bootstraps (specify 'auto' or number)"
    )
parser.add_argument(
    '--bootstrap-file',
    nargs=1,
    type=str,
    default=[None],
    dest='bootstrap_file',
    help="specify the name of a file to which all bootstrap trees should be saved"
    )
parser.add_argument(
    '-m', '--model',
    nargs=1,
    type=str,
    default=[None],
    dest='model',
    help="choose the model of sequence evolution, GTRGAMMA (default with less than 50 taxa) or GTRCAT (default with 50 taxa or more)"
    )
parser.add_argument(
    '-s', '--seed',
    nargs=1,
    type=int,
    default=[-1],
    dest='seed',
    help="random number seed"
    )
parser.add_argument(
    'infile',
    nargs='?',
    type=argparse.FileType('r'),
    default='-',
    help='the input file name'
    )
parser.add_argument(
    'outfile',
    nargs='?',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='the output file name'
    )
parser.add_argument(
    '--verbose',
    action='store_true',
    dest='verbose',
    help='write verbose output'
    )
args = parser.parse_args()
nexus = args.nexus
parsimony = args.parsimony
one_category = args.one_category
window_start_pos = args.start[0]-1
window_end_pos = args.end[0]
exclude_taxon = args.exclude_taxon[0]
partitions = args.partitions[0]
constraint = args.constraint[0]
bootstraps = args.bootstraps[0]
bootstrap_file = args.bootstrap_file[0]
outgroup = args.outgroup[0]
model = args.model[0]
seed = args.seed[0]
verbose = args.verbose
infile = args.infile
outfile = args.outfile
if infile.isatty():
    print('No input file specified, and no input piped through stdin!')
    sys.exit(0)
instring = infile.read()
inlines = instring.split('\n')

# Make sure sensible values are specified for the window start and end.
if window_start_pos < 0:
    print("ERROR: The start position of the analysis window must be at least 1!")
    sys.exit(1)
elif window_end_pos != -1:
    if window_end_pos <= window_start_pos:
        print("ERROR: The end position of the analysis window must be greater than the start position!")
        sys.exit(1)

# Make sure a sensible value is specified for the model.
if model != None:
    if model != "GTRCAT":
        if model != "GTRGAMMA":
            print("ERROR: Unknown model of sequence evolution specified: " + model + "!")
            sys.exit(1)

# Make sure the specification of bootstrap replicates is ok.
bootstraps_on = False
bootstraps_auto = False
bootstraps_number = 0
if bootstraps != None:
    bootstraps_on = True
    if bootstraps == "auto":
        bootstraps_auto = True
    elif int(bootstraps) > 0:
        bootstraps_number = int(bootstraps)
    else:
        print("ERROR: Bootstrap number could not be read!")
        sys.exit(1)

# Determine the input format, which could be phylip or fasta.
input_format = "phylip"
if inlines[0][0] == ">":
    input_format = "fasta"
elif inlines[0].lower().strip() == "#nexus":
    input_format = "nexus"
elif "[[Samples]]" in instring:
    input_format = "arlequin"

# Before writing it to temporary file, it is trimmed with -f and -t.
record_ids = []
record_seqs = []
if input_format == 'phylip':
    for inline in inlines[1:]:
        if inline != '':
            record_id = inline.split()[0]
            record_seq = inline.split()[1]
            if window_end_pos == -1:
                record_seq = record_seq[window_start_pos:]
            else:
                record_seq = record_seq[window_start_pos:window_end_pos]
            record_ids.append(record_id)
            record_seqs.append(record_seq)
elif input_format == 'fasta':
    for inline in inlines:
        if len(inline) > 0:
            if inline[0] == ">":
                record_id = inline[1:].strip().split("[")[0]
                record_ids.append(record_id)
                record_seqs.append("")
            elif inline != '':
                record_seqs[-1] += inline.strip()
elif input_format == 'arlequin':
    in_alignment = False
    for inline in inlines:
        if len(inline) > 0:
            if inline.strip() == "SampleData= {":
                in_alignment = True
            elif inline.strip() == "}":
                in_alignment = False
            elif in_alignment == True:
                record_id = inline.split()[0]
                record_seq = inline.split()[2]
                if window_end_pos == -1:
                    record_seq = record_seq[window_start_pos:]
                else:
                    record_seq = record_seq[window_start_pos:window_end_pos]
                record_ids.append(record_id)
                record_seqs.append(record_seq)
elif input_format == 'nexus':
    in_alignment = False
    for inline in inlines:
        if len(inline) > 0:
            if inline.lower().strip() == "matrix":
                in_alignment = True
            elif inline.strip() == ";":
                in_alignment = False
            elif in_alignment == True:
                record_id = inline.split()[0]
                record_seq = inline.split()[1]
                if window_end_pos == -1:
                    record_seq = record_seq[window_start_pos:]
                else:
                    record_seq = record_seq[window_start_pos:window_end_pos]
                record_ids.append(record_id)
                record_seqs.append(record_seq)
else:
    print("ERROR: Unknown input format:", input_format, "!")
    sys.exit(1)

# If an outgroup has been specified, make sure it's found among the record ids.
if outgroup != None:
    outgroup_found = False
    for record_id in record_ids:
        if outgroup == record_id:
            outgroup_found = True
    if outgroup_found == False:
        # See whether one of the record ids includes the outgroup name.
        outgroup_match_found = False
        for record_id in record_ids:
            if outgroup in record_id:
                outgroup_match_found = True
                outgroup = record_id
                print("WARNING: outgroup id changed to " + record_id + " as " + outgroup + " could not be found.")
        if outgroup_match_found == False:
            print("ERROR: The outgroup could not be found.")
            sys.exit(1)

# If -x was specified, remove all records for which the name matches the specified string.
if exclude_taxon != None:
    tmp_record_ids = []
    tmp_record_seqs = []
    for x in range(len(record_ids)):
        if exclude_taxon not in record_ids[x]:
            tmp_record_ids.append(record_ids[x])
            tmp_record_seqs.append(record_seqs[x])
    record_ids = tmp_record_ids
    record_seqs = tmp_record_seqs

# Make sure all record_seqs are of the same length.
for record_seq in record_seqs[1:]:
    if len(record_seq) != len(record_seqs[0]):
        print('WARNING: Not all sequences are of the same length!')
# Get the maximum length of record ids.
max_record_id_length = 0
for record_id in record_ids:
    if len(record_id) > max_record_id_length:
        max_record_id_length = len(record_id)
# Test whether all sequences are uninformative.
# Write the phylip string.
phylip_string = str(len(record_ids)) + ' ' + str(len(record_seqs[0])) + '\n'
for x in range(len(record_ids)):
    phylip_string += record_ids[x].ljust(max_record_id_length + 2) + record_seqs[x] + '\n'

# Set parameters.
if seed == -1:
    seed1 = str(random.randint(1, 100000))
else:
    seed1 = str(seed)
    random.seed(seed)
seed2 = str(random.randint(1, 100000))
run_id = "run_" + str(random.randint(100000, 200000))
if model == None:
    if len(record_ids) < 50:
        model = "GTRGAMMA"
    else:
        model = "GTRCAT"
if one_category == True:
    model = "GTRCAT"

# Write the alignment to a temporary file, so that RAxML can run even if data is piped into this script.
tmp_in = tempfile.NamedTemporaryFile(delete=False)
tmp_in.write(phylip_string.encode('utf-8'))
tmp_in.close()

# Run raxml.
FNULL = open(os.devnull, 'w')
call_list = ["raxml", "-T",  "6",  "-s", tmp_in.name, "-n", run_id, "-m", model, "-p", seed1, "-O", "--no-bfgs"]
if bootstraps_on == True:
    call_list.append("-f")
    call_list.append("a")
    call_list.append("-x")
    call_list.append(seed2)
    call_list.append("-N")
    if bootstraps_auto == True:
        call_list.append("autoMRE")
    else:
        call_list.append(str(bootstraps_number))
if outgroup != None:
    call_list.append("-o")
    call_list.append(outgroup)
if parsimony == True:
    call_list.append("-y")
if partitions != None:
    call_list.append("-q")
    call_list.append(partitions)
if constraint != None:
    call_list.append("-g")
    call_list.append(constraint)
if one_category == True:
    call_list.append("-V")
if verbose:
    call(call_list)
else:
    call(call_list, stdout=FNULL, stderr=FNULL)

# Read the ML tree file.
if bootstraps_on == True:
    if os.path.isfile("RAxML_bipartitions." + run_id):
        tmp_out = open("RAxML_bipartitions." + run_id)
        tree_string = tmp_out.read()
        tmp_out.close()
    elif parsimony == True and os.path.isfile("RAxML_parsimonyTree." + run_id):
        tmp_out = open("RAxML_parsimonyTree." + run_id)
        tree_string = tmp_out.read()
        tmp_out.close()        
    else:
        print('ERROR: No RAxML tree file found!')
        sys.exit(1)
else:
    if os.path.isfile("RAxML_bestTree." + run_id):
        tmp_out = open("RAxML_bestTree." + run_id)
        tree_string = tmp_out.read()
        tmp_out.close()
    elif parsimony == True and os.path.isfile("RAxML_parsimonyTree." + run_id):
        tmp_out = open("RAxML_parsimonyTree." + run_id)
        tree_string = tmp_out.read()
        tmp_out.close()        
    else:
        print('ERROR: No RAxML tree file found!')
        sys.exit(1)

# Remove unneded files created by RAxML.
if os.path.isfile("RAxML_bipartitionsBranchLabels." + run_id):
    os.remove("RAxML_bipartitionsBranchLabels." + run_id)
if os.path.isfile("RAxML_bipartitions." + run_id):
    os.remove("RAxML_bipartitions." + run_id)
if os.path.isfile("RAxML_bootstrap." + run_id):
    if bootstrap_file == None:
        os.remove("RAxML_bootstrap." + run_id)
    else:
        os.rename("RAxML_bootstrap." + run_id,bootstrap_file)
if os.path.isfile("RAxML_bestTree." + run_id):
    os.remove("RAxML_bestTree." + run_id)
if os.path.isfile("RAxML_info." + run_id):
    os.remove("RAxML_info." + run_id)
if os.path.isfile("RAxML_log." + run_id):
    os.remove("RAxML_log." + run_id)
if os.path.isfile("RAxML_parsimonyTree." + run_id):
    os.remove("RAxML_parsimonyTree." + run_id)
if os.path.isfile("RAxML_result." + run_id):
    os.remove("RAxML_result." + run_id)

if nexus:
    # Reread the alignment from the temporary file created above. This is the easiest way it can be transformed to an alignment object.
    align = AlignIO.read(open(tmp_in.name,"r"),"phylip-relaxed")
    nexus_lines = ["#NEXUS"]
    nexus_lines.append("begin data;")
    nexus_lines.append("\tdimensions ntax=" + str(len(align)) + " nchar=" + str(align.get_alignment_length()) + ";")
    nexus_lines.append("\tformat datatype=dna missing=? gap=-;")
    nexus_lines.append("matrix")
    max_id_length = 0
    for record in align:
        if len(record.id) > max_id_length:
            max_id_length = len(record.id)
    for record in align:
        nexus_lines.append(record.id.ljust(max_id_length+2) + str(record.seq))
    nexus_lines.append(";")
    nexus_lines.append("end;")
    nexus_lines.append("")
    nexus_lines.append("begin trees;")
    nexus_lines.append("\ttree ml = [&U]" + tree_string.strip())
    nexus_lines.append("end;")
    nexus_lines.append("")
    nexus_string = "\n".join(nexus_lines)

    # Write the output string to file or STDOUT.
    outfile.write(nexus_string)

else:
    # Write the output string to file or STDOUT.
    outfile.write(tree_string)

