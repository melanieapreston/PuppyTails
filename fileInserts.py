import sys
import HTSeq
import re
import string
import glob
import os
import time
import difflib
from optparse import OptionParser
import argparse
import scipy.stats

class fileInserts():
    def __init__(self, filename):
        self.filename = filename
        self.number_of_tails_by_length = {}
        # self.number_of_tails_by_length[tailLength] = number of tails * heptamers at given length
        self.tailBases = {}
        # self.tailBases[ tailLength]['A'] = sum of As * heptamers for each tail for all tails at given length
        #self.tailsByLength = {} #self.tailsByLength[ tail length ] = [ 'ATTAA', 'TTTT' ect ]
        self.total_tRNAs = 0
        self.tails = {}
        # self.tails[ tail ] = heptamer number
        self.percentBases = {}
        #dict organized by length
        
    def find_most_common_tail_at_each_length(self):
        self.most_common_at_length = dict()
        for tail_length in range(0,44):
            self.most_common_at_length[tail_length] = [("", 0), ("", 0), ("", 0)]
        tails_without_RNA_ends = dict()
        for atail in self.tails:
            tail = atail[1]
            tail = re.sub("\.", "", tail)
            tail_length = len(tail)
            if tail in tails_without_RNA_ends:
                tails_without_RNA_ends[tail] += self.tails[atail]
            else:
                tails_without_RNA_ends[tail] = self.tails[atail]
        for tail in tails_without_RNA_ends:
            tail_length = len(tail)
            min_known_at_this_length = self.most_common_at_length[tail_length][1][-1]
            if(tails_without_RNA_ends[tail] > min_known_at_this_length):
                self.most_common_at_length[tail_length].append((tail, tails_without_RNA_ends[tail]))
                self.most_common_at_length[tail_length] = sorted(
                    self.most_common_at_length[tail_length],
                    key=lambda tup: tup[1], reverse=True)
                del self.most_common_at_length[tail_length][-1]
        self.tails_without_RNA_ends = tails_without_RNA_ends
                
    def write_table_of_most_common_tails(self, control_dataset=False, get_pvalues=False):
        self.find_most_common_tail_at_each_length()
        if(control_dataset and get_pvalues):
            pvalues = self.determine_significance_of_most_common_tails(control_dataset)
            significance_by_length = self.find_significance_by_length(control_dataset)
        li = "Length"
        li += "\tTail\tNumber\tp value" * 3
        li += "\tp value of this length\n"
        for tail_length in range(0,44):
            li += "%i" % tail_length
            for a_tail_tup in self.most_common_at_length[tail_length]:
                li += "\t%s\t%.0f" % (a_tail_tup[0], a_tail_tup[1])
                if(get_pvalues and a_tail_tup[0] != ''):
                    li += "\t%e" % pvalues[a_tail_tup[0]]
                if(get_pvalues and a_tail_tup[0] == ''):
                    li += "\t1"
            if(get_pvalues):
                li += "\t%e" % significance_by_length[tail_length]
            li += "\n"
        return li
    
    def find_significance_by_length(self, control_dataset):
        significance_by_length = dict()
        total = sum(x for x in self.number_of_tails_by_length.values())
        total_control = sum(x for x in control_dataset.number_of_tails_by_length.values())
        for length in range(0,44):
            if(length in self.number_of_tails_by_length):
                n_obs = self.number_of_tails_by_length[length]
            else:
                n_obs = 0
            if(length in control_dataset.number_of_tails_by_length):
                n_obs_in_control = control_dataset.number_of_tails_by_length[length]
            else:
                n_obs_in_control = 0
            table = [
                    [n_obs, total-n_obs],
                    [n_obs_in_control, total_control - n_obs_in_control]]
            oddsratio, pvalue = scipy.stats.fisher_exact(table)
            significance_by_length[length] = pvalue
            if(n_obs + n_obs_in_control < 10):
                pvalue = 1.0
        return significance_by_length
    
    def determine_significance_of_most_common_tails(self, control_dataset):
        print "%s vs control %s" % (self.filename, control_dataset.filename)
        try:
            control_dataset.tails_without_RNA_ends[0]
        except:
            control_dataset.find_most_common_tail_at_each_length()
        pvalues = dict()
        total = sum(x for x in self.tails_without_RNA_ends.values())
        total_control = sum(x for x in control_dataset.tails_without_RNA_ends.values())
        for tail_length in self.most_common_at_length:
            for a_tup in self.most_common_at_length[tail_length]:
                tail = a_tup[0]
                n_obs = a_tup[1]
                if(tail == ''):
                    continue
                #total = float(self.total_tRNAs)
                if tail in control_dataset.tails_without_RNA_ends:
                    n_obs_in_control = control_dataset.tails_without_RNA_ends[tail]
                else:
                    n_obs_in_control = 0
                #total_control = control_dataset.total_tRNAs
                table = [
                    [n_obs, total-n_obs],
                    [n_obs_in_control, total_control - n_obs_in_control]]
                oddsratio, pvalue = scipy.stats.fisher_exact(table)
                pvalues[tail] = pvalue
        return pvalues
                
    def add_args(self, args):
        self.onlyCCA = args.onlyCCA
        self.noNegativeEnrichment = args.noNegativeEnrichment
        self.percents = args.percents
        
    def normalize(self):
        """Only normalizes the number_of_tails_by_length dict.
        """
        self.total_tRNAs = float(sum(self.number_of_tails_by_length.values()))
        if(self.total_tRNAs == 0):
            print "Error! 0 tRNAs in %s" % (self.filename)
            return False        
        for length in self.number_of_tails_by_length:
            self.number_of_tails_by_length[length] = float('1e6') * float(self.number_of_tails_by_length[length])/self.total_tRNAs
        print "Total tails = %f" % float(self.total_tRNAs)
        print "After normalization = %f" % float(sum(self.number_of_tails_by_length.values()))

    def subtractControl(self, control, src_path):
        for t in control.tails:
            if t in self.tails:
                self.tails[t] = self.tails[t] - control.tails[t]            
        self.filename =  src_path + "/controlSubtracted/%s" % (os.path.basename(self.filename))
        subF = open(self.filename, "w")
        for t in self.tails:
            st = t.split("|")
            line = "%s\t%s\t%i\n" % (st[0], st[1], self.tails[t])
            subF.write(line)
        subF.close()
        
    def add_tail(self, tRNAend, tail, numHeptamersWithTail):
        tRNA_tail = (tRNAend, tail)
        if(self.onlyCCA):
            if(tRNAend != "CCA"):
                return
        if(tRNA_tail in self.tails):
            self.tails[tRNA_tail] += numHeptamersWithTail
        else:
            self.tails[tRNA_tail] = numHeptamersWithTail
            
    def process_tails(self, num_samples_averaged_together):
        for t in self.tails:
            #self.tails[t] = float(self.tails[t])/float(num_samples_averaged_together)
            #self.tails[t] = int(self.tails[t])
            numHeptamersWithTail = self.tails[t]
            tail = t[1]
            tRNAend = t[0]
            tail = re.sub("\.", "", tail)
            if(self.noNegativeEnrichment):
                numHeptamersWithTail = max(numHeptamersWithTail, 0)
            # Add to dict of tails by length counting heptamers.
            if len(tail) in self.number_of_tails_by_length:
                self.number_of_tails_by_length[len(tail)] += numHeptamersWithTail
            else:
                self.number_of_tails_by_length[len(tail)] = numHeptamersWithTail            
            # Count the bases.
            if len(tail) in self.tailBases:
                s = 0
                for b in self.tailBases[len(tail)]:
                    s += self.tailBases[len(tail)][b]     
                self.tailBases[len(tail)] = {
                    'A': self.tailBases[len(tail)]['A'] + (numHeptamersWithTail * tail.count('A')),
                    'C': self.tailBases[len(tail)]['C'] + (numHeptamersWithTail * tail.count('C')),
                    'G': self.tailBases[len(tail)]['G'] + (numHeptamersWithTail * tail.count('G')),
                    'T': self.tailBases[len(tail)]['T'] + (numHeptamersWithTail * tail.count('T')),
                    'N': self.tailBases[len(tail)]['N'] + (numHeptamersWithTail * tail.count('N'))}
            else:
                self.tailBases[len(tail)] = {
                    'A': (numHeptamersWithTail * tail.count('A')),
                    'C': (numHeptamersWithTail * tail.count('C')),
                    'G': (numHeptamersWithTail * tail.count('G')),
                    'T': (numHeptamersWithTail * tail.count('T')),
                    'N': (numHeptamersWithTail * tail.count('N'))}
            # Has adding this tail caused an arithmetic error yet?
            s=0
            for b in self.tailBases[len(tail)]:
                s += self.tailBases[len(tail)][b]
            
            if((self.number_of_tails_by_length[len(tail)] * len(tail)) == s):
                pass
            else:
                print "number_of_tails_by_length array holds %i " % (self.number_of_tails_by_length[ len(tail) ])
                es = (self.number_of_tails_by_length[ len(tail) ] * len(tail))
                #print "exepct base sum to be %i, base sum actually %i. sequence %s" % (es, s, tail)
                #print "did not match (heptamer number %i at tail length %i * tail length) == sum of bases %i)" % (numHeptamersWithTail, len(tail), s)
                #raw_input("Press Enter to continue...")

    def calculate_bases(self):        
        for tailLen in self.number_of_tails_by_length:
            if (tailLen == 0):
                continue
            numBases = 0
            for base in self.tailBases[ tailLen ]:
                numBases = numBases + self.tailBases[ tailLen ][base]
            #numBases = self.tailBases[ tailLen ]['A'] + self.tailBases[ tailLen ]['C'] + self.tailBases[ tailLen ]['G'] + self.tailBases[ tailLen ]['T'] + self.tailBases[ tailLen ]['N']
            if(numBases):
                self.percentBases[ tailLen ] = {}
                for base in self.tailBases[ tailLen ]:
                    self.percentBases[ tailLen ][ base ] = float(self.tailBases[ tailLen ][base])/float(numBases)

    def write_tails_file(self, outfname):
        tF = open(outfname, "w")
        for t in self.tails:
            outLine = "%s\t%s\t%i\n" % (t[0], t[1], self.tails[t])
            tF.write(outLine)
        outLine = r'sort -k1,2 ' + outfname + ' > sorted.tmp'
        print outLine
        os.system(outLine)
        outLine = r'mv sorted.tmp ' + outfname
        print outLine
        os.system(outLine)

    def write_for_graph(self, graF):
#       graF.write("Data for file %s\n" % (self.filename))
        #output the number of tails at given length
        graF.write("number tails in %s:\t" % (self.filename))
        for aLen in range(0, 44):
            if(aLen in self.number_of_tails_by_length):
                outLine = "%i\t" % self.number_of_tails_by_length[aLen]
                graF.write(outLine)
            else:
                graF.write("0\t")
        graF.write("\n")
        if(self.percents):
        #output the number of bases in tails at the given length
            for base in ['A', 'G', 'C', 'T', 'N']:
                graF.write("%s:\t" % (base) )
                for aLen in range(0, 44):
                    if(aLen in self.percentBases):
                        outLine = "%f\t" % self.percentBases[aLen][base]
                        graF.write(outLine)
                    else:
                        graF.write("0\t")
                graF.write("\n")
        else:            
            #output the number of bases in tails at the given length
            for base in ['A', 'G', 'C', 'T', 'N']:
                graF.write("%s:\t" % (base) )
                for aLen in range(0, 44):
                    if(aLen in self.tailBases):
                        outLine = "%i\t" % self.tailBases[ aLen ][ base ]
                        graF.write(outLine)
                    else:
                        graF.write("0\t")
                graF.write("\n")
        #output the axis for the graph
        graF.write("tail length:\t")
        for aLen in range(0,44):
            graF.write("%i\t" % (aLen))
        graF.write("\n")
        #output the sequence diversity for the graph: the number of unique sequences seen
        #calculate the sequence diversity first:
        uniqueSeqsByLen = dict()
        for fulltail in self.tails:
            t = fulltail[1]
            if(len(t) in uniqueSeqsByLen):
                uniqueSeqsByLen[ len(t) ].add(t)
            else:
                uniqueSeqsByLen[ len(t) ] = set( [t] )
        if(1 in uniqueSeqsByLen):
            uniqueSeqsByLen[ 1 ] = uniqueSeqsByLen[ 1 ] - set(['.']) #get rid of the '.' tails from tail length one
        #then output the sequence diversity:
        graF.write("Number of unique tail sequences\t")
        for aLen in range(0, 44):
            if(aLen in uniqueSeqsByLen):
                outLine = "%i\t" % (len(uniqueSeqsByLen[ aLen ]))
                graF.write(outLine)
            else:
                graF.write("0\t")
        graF.write("\n")
        #output the axis for the graph (again)
        graF.write("tail length:\t")
        for aLen in range(0,44):
            graF.write("%i\t" % (aLen) )
        graF.write("\n")
