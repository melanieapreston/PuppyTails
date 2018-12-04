import sys
import HTSeq
import re
import string
import glob
import os
import time
import itertools
import difflib
from optparse import OptionParser
import timeit
import argparse

"""
This program has been modifed to find the end sequence of the Y3H RNA, tRNA,
or any given mRNA. If not using the default tRNA end, the
targetted RNA end is given in a file passed with the --seq_file <filename> argument.
For a tRNA, we expect, with X = <some sequence>,
X2 (optional) = <some other sequence>, ect:
file format:
X
XC
XCC
XCCA
X2C
X2CC
X2CCA
ect.
The X is the fixed sequence. The variations after the fixed sequence are
the variable sequence.
This program will open every .fastq file in the current directory.
It will extract tRNA 3'end sequences and generate a .inserts file for
each .fastq file, containing the tRNA end, tail and random heptamer for
each tRNA end.
It will then read each .inserts file and remove duplicates, using the
random heptamer. It generates a .inserts.noDups file that contains the
tRNA end, tail and number of tRNA ends (with different heptamers) containing
that end type. The .inserts.noDups file is sorted by tRNA end, and then by
length of insert.
 i.e.
tRNA end \t tail \t number of times this RNA is found with a unique heptamer
A \t TTT \t 3
"""


def parse_args():
    """
    Determine if tRNAs with no tails are included (DEFAULT: included).
    """
    parser = argparse.ArgumentParser(description=
    """Meta-analysis of CLIP-seq .peaks files.""")
    parser.add_argument('-b', '--blanks', default=True, action='store_false',
                        help="""(Optional)  tRNAs with no tails are not counted.""")
    parser.add_argument('-s', '--seq_file', default=False,
                        help=""" File of tRNA tails to search through.
If not given, will default to searching for CGACAAC.""")
    parser.add_argument('-d', '--directory', default="./",
                        help="""Directory of fastq files""")
    parser.add_argument('-y', '--y3h_mode', default=False,
                        help="""Tails are Y3H tails.""")
    args=parser.parse_args()
    return args


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)


def rc(s): 
    s = s[::-1] 
    s = complement(s) 
    return s


class RNAend():
    def __init__(self, heptamer, variable_seq, insert):
        self.heptamer = heptamer
        if(variable_seq):
            self.variable = variable_seq
        else:
            self.variable = "."
        if(insert):
            self.insert = insert
        else:
            self.insert = "."
        self.count = 1
        self.checkedForDups = 0
        self.isAduplicate = 0
    def objects(self):
        return "\t".join([self.variable, self.insert, self.heptamer])
    def write(self):
        print self.heptamer, "\t", self.variable,"\t",self.insert, "\t", self.count


def trim_RNA_end_or_identify_other_RNA(R1seq, fixed_seq):
    mismatches = [x == y for (x,y) in zip(R1seq[0:41], fixed_seq)].count(False)
    #print "Found %d mismatches between\n%s and \n%s" % (mismatches, R1seq[0:41], tRNA_seq)
    if(mismatches>2):
        return False #This is not a tRNA
    else:
    #if this is a tRNA, get the bit after the tRNA end
        #print "end = %s" % (R1seq[42:])
        return R1seq[42:]


def is_this_adapter(R2seq):
    rcSeq = rc(R2seq)
    adapterPats = {}
    adapterPats['pat1'] = re.search(r'TGGAATTCT', rcSeq)
    for key in adapterPats:
        if adapterPats[key] is not None:
            return True # see adapter sequence
        else:
            pass
    return False #don't see adapter


def parse_input_RNA_ends_file(fname):
    if not fname:
        shortest_seq = rc('CGACAAC')
        variable_seqs = ['CCA', "CC", "C", "."]
        #variable_seqs = ['TGG', "GG", "G", "."]
        return (shortest_seq, variable_seqs)
    seqs_from_file = list()
    with open(fname, 'r') as f:
        for li in f:
            if(re.match('\A\s*$', li)):
                continue
            li = li.rstrip('\n')
            if(re.match('[^ACGTacgtuU]', li)):
                print "Unexpected sequence format in %s" % li
                continue
            seqs_from_file.append(li.rstrip('\n'))
    fixed_seq = ""
    variable_seqs = []
    shortest_seq = seqs_from_file[0]
    for aseq in seqs_from_file:
        if len(aseq) < len(shortest_seq):
            shortest_seq = aseq
    print "seqs: %s" % str(seqs_from_file)
    print "shortest: %s" % shortest_seq
    for aseq in seqs_from_file:
        mismatches = [x == y for (x,y) in zip(aseq, shortest_seq)].count(False)
        if(mismatches > 0):
            print """
Error: all sequences in the input RNA ends file should begin with and
entirely include the shortest sequence. That is, the expected format is:
AGAGA
AGAGAC
AGAGACC
ect.
In this case, sequence %s does not include sequence %s.
""" % (aseq, shortest_seq)
            sys.exit()
        if(len(aseq) > len(shortest_seq)):
            variable_seqs.append(aseq[len(shortest_seq):])
        else:
            variable_seqs.append('.')
    print variable_seqs
    return (shortest_seq, variable_seqs)


def check_if_new_R1_file_get_R2_filename(filename, R1_files_parsed, R2_files_parsed):
    filenamePat = r'(\w*\.*[\w|-]*)_R(\d).fastq'
    fastq_dir = os.path.dirname(filename)
    fastq_basename = os.path.basename(filename)
    pat = re.search(filenamePat, fastq_basename)
    if(pat.group(2) == "1"):
        R2_filename = fastq_dir + '/' + pat.group(1) + "_R2.fastq"
        # This is R1. Expect R2 to be file <R2_filename>
        if(filename in R1_files_parsed):
	# We have parsed this file before.
            return False
        else:
            if(R2_filename in R2_files_parsed):
                return False
            else:
                return R2_filename
    if(pat.group(2) == "2"):
	# Only parse R2 files with R1 - skip this file for now
        return False


def parse_R1_file(R1_fastq_file, fixed_seq):
    R1_reads_with_fixed_seq = list()
    R1_reads_trimmed = list()
    cat = {
        'numRead': 0,
        'numRep': 0,
        'numNoTail': 0,
        'numHasTail': 0}
    for aRead in R1_fastq_file:
        cat['numRead'] += 1
        readStr = aRead.seq
        # First check if this read is a tRNA. If not, discard. If it is, keep
        # the sequence downstream of the tRNA end.
        post_fixed_seq = trim_RNA_end_or_identify_other_RNA(aRead.seq, fixed_seq)
        if(post_fixed_seq):
            R1_reads_with_fixed_seq.append(cat['numRead'])
            R1_reads_trimmed.append(post_fixed_seq)
        else:
            continue # Not a tRNA. Go to next read.
    fraction_fixed = float(len(R1_reads_with_fixed_seq))/float(cat['numRead'])
    outLine = """From R1 file:\nReads:\t%d\tWith RNA end:%d\tFraction:%f\n**""" % (
        cat['numRead'], len(R1_reads_with_fixed_seq),
                    fraction_fixed)
    print outLine
    return (R1_reads_with_fixed_seq, R1_reads_trimmed, outLine)


def is_this_pre_tRNA(R2seq):
    # Heptamer has not been cut off.
    rcSeq = rc(R2seq) #make it easier to think, using the reverse complement
    pre_tRNA_pats = {}
    pre_tRNA_pats['pat_upstream'] = re.search(r'TTTCCGAAAT', rcSeq)
    pre_tRNA_pats['pat_downstream'] = re.search(r'ATCATTTTTT', rcSeq)
    pre_tRNA_pats['pat_downstream2'] = re.search(r'TTGTTGCGAGTTGT', rcSeq)
    pre_tRNA_pats['pat3'] = re.search(r'CCGAAATTTCG', rcSeq)
    for key in pre_tRNA_pats:
        if pre_tRNA_pats[key] is not None:
            return True #this is a pre-tRNA sequence
    return False #this is not a pre-tRNA sequence


def parse_R2_file(R2_fastq_file, outF, counter, fixed_seq, variables,
                  R1_reads_with_fixed_seq):
    numRead = 0
    variable_ends = variables
    #fixed_end_seq_pat = r'\A(\w+)GAAAGAT'
    seq_pats = list()
    for a_variable in variables:
        if(a_variable == '.'):
           a_variable = ""
        seq_pats.append(r'\A(\w+)(' + rc(a_variable) + ')(' + rc(fixed_seq) + ')')
    for numRead, R2_read in enumerate(R2_fastq_file, start=1):
        counter['num_reads'] += 1
        prevCounter = counter.copy()
        if(is_this_adapter(R2_read.seq)):
            counter['adapter'] += 1
            continue
        if(is_this_pre_tRNA(R2_read.seq)):
            counter['pre-tRNA'] += 1
            continue
        counter['not_adapter_or_pre_tRNA'] += 1
        found_end = False
        for a_pat in seq_pats:
            if found_end:
                continue
            end_match = re.search(a_pat, R2_read.seq)
            if ((end_match is not None) and (len(end_match.group(1)) >= 7)):
                found_end = True
                 # Skip sequences with heptamers < 7
                counter['R2_tRNA'] += 1
                heptamer = end_match.group(1)[0:7]
                found_var = rc(end_match.group(2))
                found_fixed = end_match.group(3)
                tail = end_match.group(1)[7:]
                if(tail):
                    tail = rc(tail)
                    counter['has_tail'] += 1
                    newEnd = RNAend(heptamer, found_var, tail)
                    outF.write(newEnd.objects())                        
                    outF.write("\n")
                else:
                    counter['no_tail'] += 1
                    #print "no_tail counter: %i" % counter['no_tail']
                    if(includeBlankTails):
                        newEnd = RNAend(heptamer, found_var, '')
                        outF.write(newEnd.objects())
                        outF.write("\n")
        if(not found_end):  # Don't see a tRNA end.
            try:
                indx = R1_reads_with_fixed_seq.index(numRead)
                #most of these look like imperfect matches to the R2 regex pattern.
                #so we try an alternative:
                newEnd = RNAend(R2_read.seq[0:7], '?', rc(R2_read.seq[7:]))
                outF.write(newEnd.objects())
                outF.write("\n")
                counter[ 'R1_tRNA_only' ] += 1
                #print R2_read.seq[7:]
            except:
                counter['no_tRNA_in_either_read'] += 1


def remove_duplicates(inserts_filename, no_dups_filename):
        """Remove duplicates by assuming the file is sorted
        by tRNA end, and then by tail
        #It will fail if the file is not sorted.
        #A "hand" is defined comprising the most recent tRNA end and tail, and the number of
        #reads with distinct random barcodes is defined as the size of the hand.
        #When a line (AKA read) is hit with a different tRNA end or tail, the hand is written
        #to the output file and a new hand is defined based on the current line.
        #Currently any difference between heptamers counts as a unique sequence.
        #As a result, unique reads are over-estimated proportional to the number
        #of times the same PCR is sequenced (that is, sequencing error will result
        #in more 'unique' reads by misreading the same heptamer).
        #4^7 = 16383 possible unique heptamers
        """
        cmd = "sort -s -k1,2 %s > %s_sorted" % (inserts_filename, inserts_filename) 
        os.system(cmd)
        cmd = "mv %s_sorted %s" % (inserts_filename, inserts_filename) 
        os.system(cmd)
        print """Removing duplicates...
This can be very slow for files with a great deal of duplication.
For speed, reads with blank tails can be excluded by running with -b."""
        insertsF = open(inserts_filename, 'r')
        noDupsF = open(no_dups_filename, 'w')
        tRNAinHand = "."
        tailInHand = "."
        hand = set()
        first_seq = True
        start = time.time()
        for li in insertsF:
            #print "cf " + l + " with hand -" + tRNAinHand + "-\t-" + tailInHand +"-"
            #lines are tRNA\tTail\tHeptamer
            s = li.rstrip("\n").split('\t')
            if((s[0] == tRNAinHand) and (s[1] == tailInHand)):
                hand.add(s[2])
                #print "Adding to hand %s. Hand is now size=%i" % (li, len(list(hand)))
                #print "Hand is %s variable and %s tail" % (tRNAinHand, tailInHand)
            else:
                if(not first_seq):
                    num_heptamers = str(len(list(hand)))
                    outLine = "\t".join([tRNAinHand, tailInHand, num_heptamers]) + "\n"
                    #print "Printing to no dups... %s" % outLine
                    noDupsF.write(outLine)
                first_seq = False
                hand = set()
                hand.add(s[2])
                tRNAinHand = s[0]
                tailInHand = s[1]
        insertsF.close()
        noDupsF.close()
        end = time.time()
        timeTaken = end - start
        print "Time taken to remove duplicates:%f" % (timeTaken)
        #sort the no-dups inserts file by tRNA end and then by insert length
        #this bit here will not work if there is a huge number of different inserts, beyond memory capacity
        #this can be replaced later by a different method if there turns out to be a memory problem
        print "Sortin the no-duplicates file..."
        sortedData = sorted(open(no_dups_filename, 'r').readlines(), key = lambda line: len(line.split("\t")[1]) )
        sortedData = sorted(sortedData , key = lambda line: line.split("\t")[0] )
        open(no_dups_filename, 'w').writelines(sortedData)
        summaryOutFile.close()
        cmd = "sort -s -k3nr %s > %s" % (no_dups_filename, inserts_filename)
        print cmd
        os.system(cmd)
        os.system("rm %s" % no_dups_filename)
        #os.system("rm %s_sorted" % no_dups_filename)

           
if __name__ == '__main__':
    src_path = os.path.dirname(os.path.realpath(__file__))
    args = parse_args()
    includeBlankTails = args.blanks
    if(includeBlankTails):
        print "RNAs with no tail are included"
    else:
        print "RNAs with no tail are excluded"
    (fixed_seq, variables) = parse_input_RNA_ends_file(args.seq_file)
    print """Fixed seq end=%s\nVariables=%s""" % (fixed_seq, str(variables))
    R1_files_parsed = []
    R2_files_parsed = []
    #wipe any existing results summary file
    summaryOutFile = open('results.txt', mode='w')
    summaryOutFile.close()
    # This is the end of the full-length Y3H RNA sequence
    #tRNA_seq = r'GATCGGAATTCCCCCATATCCAACTTCCAATTTAATCTTTCTTTT'
    #R1_fixed_seq = r'GATCGGAATTCCCCCATATCCAACTTCCAATTTAATCTTTCTTTT'
    R1_fixed_seq = r'GAGGATCACCCATGTCGCAGGTTCGAGTCCTGCAGTTGTCG'
    if args.y3h_mode:
        R1_fixed_seq = r'GTctgcaggtcgactctagaAAACATGAGGA'.upper()
    for filename in glob.glob(args.directory.rstrip('/') + '/*.fastq'):
        bytRNAfrag = {}
        R2_filename = check_if_new_R1_file_get_R2_filename(
            filename, R1_files_parsed, R2_files_parsed)
        if(not R2_filename):
            continue
        R2_fastq_file = HTSeq.FastqReader(R2_filename, 'solexa')
        li = """\n***\nR1:%s, R2:%s\n""" % (filename, R2_filename)
        summaryOutFile = open('results.txt', 'a')
        summaryOutFile.write(li)
        print li.rstrip('\n')
        base_filename = os.path.basename(filename)
        if(not os.path.exists(src_path + '/inserts')):
            os.system('mkdir ' + src_path + '/inserts')
        inserts_filename = src_path + '/inserts/' + base_filename + ".inserts"
        outF = open(inserts_filename, mode='w')
        #open the R1 fastq files
        R1_fastq_file = HTSeq.FastqReader(filename,'solexa')
        (R1_reads_with_fixed_seq, R1_reads_trimmed, outLine) = parse_R1_file(
            R1_fastq_file, R1_fixed_seq)
        summaryOutFile.write(outLine)
        counter = {'R2_tRNA': 0, 'R1_tRNA_only': 0,
                    'no_tRNA_in_either_read': 0, 'pre-tRNA': 0,
                    'adapter': 0, 'no_tail': 0,
                    'has_tail': 0, 'num_reads': 0,
                    'not_adapter_or_pre_tRNA': 0}
        parse_R2_file(R2_fastq_file, outF, counter, fixed_seq, variables,
                      R1_reads_with_fixed_seq)
        fraction = {}
        for key in sorted(counter.iterkeys(), key=lambda k: counter[k], reverse=True):
            fraction[key] = float(counter[key])/float(counter['num_reads'])
            outLine = "%s:\t%d\tFraction of reads:%f\n" % (key, counter[key], fraction[key])
            summaryOutFile.write(outLine)
            print outLine.rstrip('\n')
        #fraction_see_tRNA = float(num_has_tRNA)/float(numRead)
        #fraction_no_tRNA_in_either = float(no_tRNA_in_either)/float(numRead)      
        outF.close()
        #now sort
        no_dups_filename = inserts_filename + '.noDups'
        remove_duplicates(inserts_filename, no_dups_filename)


    #rcThree = "CCTTGGCACCCGAGAATTCCA"
    #rctRNA = "GAAAGATTAAATTGGAAGTTGGATATGGGGGAATTCCGATC"
    #tRNA_end_rc = r'\A(\w+)GAAAGATTAAATTG' 
    #r1seq = r'\A(\w*)TGGAA(\w+)' 
    #t = r'\A(\w*)CTTTC(\w+)'
    #findBothSeqs = re.compile(tRNA_end_rc)
