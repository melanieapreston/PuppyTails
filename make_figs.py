import sys
import HTSeq
import re
import string
import glob
import os
import time
import difflib
import argparse


def parse_input():
    parser = argparse.ArgumentParser(description="""
    USAGE: python make_figs.py -f data_file
    """)

    # If the -b option is used, tRNAs with no tails are not counted.
    # This speeds up the removal of duplicates for large datasets
    #parser.add_option("-b", "--blanks", action="store_false", dest="includeBlankTails", default=True)

    parser.add_argument("-f", "--data_file", action="store",
                        dest="data_file",
                        help="Filename of data.")
    args = parser.parse_args()
    return args


def write_most_common_tails(inserts, base_filename, control=False):
    for exp in inserts:
        with open("%s_%s" % (base_filename,
                             os.path.basename(exp).rstrip('.inserts').rstrip(
                                 '.fastq')),
                  'w') as f:
            if(not control):
                lines = inserts[exp].write_table_of_most_common_tails(control)
            if(control):
                lines = inserts[exp].write_table_of_most_common_tails(
                    control, get_pvalues=True)
            f.write(lines)


def parse_data_file(filename):
    data = {}
    print "Opening %s with file size %i..." % (
        filename, os.path.getsize(filename))
    with open(filename, 'r') as f:
        dataset = ""
        for li in f:
            #print li
            s = li.strip('\n').split('\t')
            m = re.match(r'number tails in ([^:]+):.*', li)
            if(m is not None):
                dataset = m.group(1)
                dataset = os.path.basename(dataset)
                cur_dataset = dataset
                data[dataset] = {'n_tails': s[1:]}
                continue
            m = re.match(r'([AGCTN]):.*', s[0])
            if(m is not None):
                data[dataset][m.group(1)] = s[1:]
                continue
            m = re.match(r'tail length:.*', li)
            if(m is not None):
                data[dataset]['tail_len'] = s[1:]
                continue
            m = re.match(r'.*Number of unique.*', li)
            if(m is not None):
                data[dataset]['n_unique'] = s[1:]
                continue
    return data
            

def check_data_agreement(data):
    for exp in data:
        max_range = min(len(data[exp]['n_tails']),
                        len(data[exp]['tail_len']),
                        len(data[exp]['n_unique']))
        n_tails = 0
        for index in range(1, max_range-1):
            try:
                n_tails += float(data[exp]['n_tails'][index])
            except:
                print "Error at %s, %i" % (exp, index)
        print "%s: total tails=%f" % (exp, n_tails)
        

def write_for_R(data, src_path):
    src_path = os.path.dirname(os.path.realpath(__file__))
    files_for_R = list()
    check_data_agreement(data)
    for exp in data:
        with open("%s/figs/%s.forR" % (
            src_path, exp.rstrip('.fastq.inserts')
            ), 'w') as f:
            li = "tail_len\tn_tails\tn_unique\tA\tC\tT\tG\n"
            max_range = min(len(data[exp]['n_tails']),
                            len(data[exp]['tail_len']),
                            len(data[exp]['n_unique']))
            for index in range(0, max_range):
                li += "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    data[exp]['tail_len'][index],
                    data[exp]['n_tails'][index],
                    data[exp]['n_unique'][index],
                    data[exp]['A'][index],
                    data[exp]['C'][index],
                    data[exp]['T'][index],
                    data[exp]['G'][index])
            f.write(li)
            files_for_R.append("%s/figs/%s.forR" % (
                src_path, exp.rstrip('.fastq.inserts')))
    return files_for_R


def r_script_for_barplot(files_for_R, src_path):
    for filename in files_for_R:
        li = """
    f = read.table("%s", head=T)""" % filename
        li += """
    bases = as.data.frame(cbind(f$A, f$C, f$T, f$G))
    m = as.matrix(bases)
    outfname = "%s/figs/barplot_%s.eps"
    """ % (src_path, os.path.basename(filename))
        li += r'''
    library(RColorBrewer)
    my_cols <- brewer.pal(4, "RdBu")
    setEPS(width=5,height=3); postscript(outfname)
    barplot(t(m), xlab = 'Tail length',
    ylab = 'Percent base composition',
    legend=c('A','C','T','G'), col=my_cols)
    dev.off()
    '''
        li += """
    outfname = "%s/figs/plot_%s.eps"
""" % (src_path, os.path.basename(filename))
        li += r'''
    library(RColorBrewer)
    my_cols <- brewer.pal(4, "RdBu")
    setEPS(width=5,height=10); postscript(outfname)
    par(mfrow=c(3,1))
    plot(f$n_tails, x=f$tail_len, type='l', xlab='Tail length',
    ylab='Number of tails')
    plot(f$n_unique, x=f$tail_len, type='l', xlab='Tail length',
    ylab='Number of unique tails')
    barplot(t(m), xlab = 'Tail length',
    ylab = 'Percent base composition',
    legend=c('A','C','T','G'), col=my_cols)
    dev.off()
    '''
        with open('tmp.r', 'w') as f:
            f.write(li)
        cmdl = """R CMD BATCH tmp.r"""
        os.system(cmdl)


def make_figs(data_filename, src_path):
    print "In make_figs. Processing file %s" % data_filename
    data = parse_data_file(data_filename)
    if(not os.path.exists(src_path + "/figs")):
        print "making %s/figs" % src_path
        os.system("mkdir %s/figs" % src_path)
    files_for_R = write_for_R(data, src_path)
    r_script_for_barplot(files_for_R, src_path)

    
if __name__ == '__main__':
    src_path = os.path.dirname(os.path.realpath(__file__))
    args = parse_input()
    data = parse_data_file(args.data_file)
    if(not os.path.exists(src_path + '/figs')):
        os.system('mkdir ' + src_path + '/figs')
    files_for_R = write_for_R(data)
    r_script_for_barplot(files_for_R)
