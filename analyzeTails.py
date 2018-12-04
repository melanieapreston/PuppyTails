import sys
import HTSeq
import re
import string
import glob
import os
import time
import difflib
import argparse
from fileInserts import fileInserts
import make_figs


def parse_input(src_path):
    parser = argparse.ArgumentParser(description="""
    USAGE: python analyzeTails.py [options]
    Outputs a dataForGraph.txt file.
    """)

    # If the -b option is used, tRNAs with no tails are not counted.
    # This speeds up the removal of duplicates for large datasets
    #parser.add_option("-b", "--blanks", action="store_false", dest="includeBlankTails", default=True)
    parser.add_argument("-i", "--inserts", default="%s/inserts/" % src_path,
			help="Folder of .inserts files")
    parser.add_argument("-n", "--normalize", action="store_true",
                        dest="normalize", default=False,
                        help="normalize to per 1 million tRNAs with unique heptamer")
    #parser.add_argument("-c", "--control", action="store_true",
    #                    dest="control", default=False,
    #                    help="(Optional) Use a control. (Default: False)")
    parser.add_argument("-s", "--subtract_control", action="store_true",
                        dest="subtract_control", default=False,
                        help="(Optional) Use a control. (Default: False)")
    parser.add_argument("-c", "--control_file", action="store",
                        dest="control_file",
                        help="Filename for control subtraction.")
    parser.add_argument("-p", "--percents", action="store_true",
                        dest="percents", default=False,
                        help="output base composition as a percentage.")
    parser.add_argument("-a", "--noNegativeEnrichment", action="store_true",
                        dest="noNegativeEnrichment", default=False,
                        help="Set negative enrichments as 0 (Default: False)")
    parser.add_argument("-m", "--average_files", action="store_true",
                        dest="averageFiles", default=False,
                        help="Average the negative controls. (Default: False)")
    parser.add_argument("-o", "--only_CCA", action="store_true",
                        default=False,
                        help="Only use tails with CCA ends. (Default: False)")
    args = parser.parse_args()
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


def read_in_all_data_and_normalize(globstr, normalize, args,
                                   write_for_graph=True,
                                   data_filename="dataForGraph.txt"):
    insertsObjects = dict()
    graF = open(data_filename, mode='w')
    for filename in glob.glob(globstr):
        print "Reading in data from %s..." % filename
        with open(filename, mode='r') as tailsFile:
            insertsFile = fileInserts(filename)
            insertsFile.add_args(args)
            for line in tailsFile:
                s = line.split("\t")
                if(args.only_CCA and s[0] != "CCA"):
                    continue
                if(s[0] == "?"):
                    continue
                insertsFile.add_tail(s[0], s[1], int(s[2].rstrip('\n')))
                # Pass (tail, number of heptamers with tail).
            insertsFile.process_tails(1)
            if(normalize):
                insertsFile.normalize()
            if(write_for_graph):
                insertsFile.calculate_bases()
                insertsFile.write_for_graph(graF)
            insertsObjects[filename] = insertsFile
        print "...Read in %i tails" % len(list(insertsObjects[filename].tails))
    graF.close()
    return insertsObjects


if __name__ == '__main__':
    src_path = os.path.dirname(os.path.realpath(__file__))
    if src_path == "":
        src_path = "."
    args = parse_input(src_path)
    #includeBlankTails = options.includeBlankTails
    normalize = args.normalize
    control = args.control_file
    percents = args.percents
    args.onlyCCA = False
    inserts = dict()
    if(args.averageFiles):
        li = "Will average the following files to form the negative control."
        li += " (-m option set):"
        for f in args:
            li += " %s," % str(f)
        li = li.rstrip(',')
        li += "\nDo not use -m option with -f option."
        li += " Negatives are passed as arguments without -f flag."
        print li
    insertsObjects = {}
    if(not os.path.exists(src_path + '/tables')):
        os.system("mkdir %s/tables" % src_path)
    if(args.subtract_control):
        if(not os.path.exists(src_path + '/controlSubtracted')):
            os.system("mkdir  %s/controlSubtracted/" % src_path)
        # Read in all the data and normalize, but do no other processing.
        globstr = args.inserts + '/*.inserts'
        insertsObjects = read_in_all_data_and_normalize(globstr, normalize, args,
                                                write_for_graph=False,
                                                data_filename="%s/dataForGraph.txt" % src_path)
        if(args.averageFiles):
            # If we are averaging negatives
            controlInsertsFile = args.inserts + '/negatives_averaged.fastq.inserts'
            negativesAve = fileInserts(controlInsertsFile)
            negativesAve.add_args(args)
            for f in args:
                if(f in insertsObjects):
                    print "\n\t...Adding %s to averaged negatives object..." % f
                    for t in insertsObjects[f].tails:
                        #print "hepval=%i" % insertsObjects[f].tails[t]
                        negativesAve.add_tail(t[0], t[1], insertsObjects[f].tails[t] )
                else:
                    print "ERROR! Did not find %s as a negatives file for averaging! " % f
            #processing summed tails to finish creating this object
            negativesAve.process_tails(len(args))
            negativesAve.write_tails_file(src_path + '/negatives.tails')

            print "\n...Subtracting the averaged negatives..."
            for f in insertsObjects:
                print "\n\t...Subtracting %s from %s" % (negativesAve.filename,
                                                         insertsObjects[f].filename)
                # Writes a new tails file file.
                insertsObjects[f].subtractControl(negativesAve, src_path)
        else:
            # Not averaging negatives, using -c flag.
            # Did we find control?
            if args.control_file in insertsObjects:
                controlInsertsFile = insertsObjects[args.control_file]
                print "\n...Subtracting the control %s...." % (controlInsertsFile.filename)
                for f in insertsObjects:
                    if(insertsObjects[f].filename == args.control_file):
                        continue
                    print "\n\t...Subtracting %s from %s" % (controlInsertsFile.filename,
                                                             insertsObjects[f].filename)
                    insertsObjects[f].subtractControl(controlInsertsFile, src_path)
                    # Writes a new tails file file.          
            else:
                print "error! could not find control file in processed file objects!\n"
                sys.exit()
        # Control subtracted tails files are now in *.controlSubtracted files
        # Process them like before.
        # They are already normalized, if the -n option is set.
        globstr = src_path + '/controlSubtracted/*.inserts'
        inserts = read_in_all_data_and_normalize(globstr, normalize, args,
                        data_filename=src_path + "/controlSubtracted/dataForGraph.txt")
        make_figs.make_figs(src_path + "/controlSubtracted/dataForGraph.txt", src_path)
        make_figs.write_most_common_tails(
            inserts, src_path + '/tables/common_tails',
            control=controlInsertsFile)
    else:
        # Not subtracting control        
        globstr = (args.inserts + '/*.inserts')
        inserts = read_in_all_data_and_normalize(globstr, normalize, args,
                                                 data_filename=src_path + "/dataForGraph.txt")
        make_figs.make_figs(src_path + "/dataForGraph.txt", src_path)
        if(not args.control_file):
            control_dataset = False
        else:
            control_dataset = inserts[args.control_file]
        make_figs.write_most_common_tails(
            inserts, src_path + '/tables/common_tails',
            control=control_dataset)
