#!/usr/bin/env python
'''Simple wrapper for the LANL program codonalign
http://www.hiv.lanl.gov/content/sequence/CodonAlign/codonalign.html
It only works with fasta files and depends on Biopython.
'''
from __future__ import print_function
import sys
import os.path
import logging
import logging.handlers


from Bio import SeqIO

dn = os.path.dirname(__file__)

# Make a global logging object.
calog = logging.getLogger(__name__)
# Master sequence where stop codons are counted is always the first one
master_seq = 1


def main(filename=None, frame=None, fs_comp=None):
    '''
    Runs codon align
    '''
    import shutil
    import tempfile
    from textwrap import fill

    # set logging level
    calog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './codal.log'
    hl = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                              maxBytes=100000, backupCount=5)
    f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s\
                          %(lineno)d %(message)s")
    hl.setFormatter(f)
    calog.addHandler(hl)
    calog.info(' '.join(sys.argv))
    if len(sys.argv) == 1:
        calog.info('running with all default options')

    out_dir = tempfile.mkdtemp(dir='./')
    calog.info('Writing to directory %s' % out_dir)

    # convert to table format: id <TAB> sequence
    table_file = '%s/table.txt' % out_dir
    oh = open(table_file, 'w')
    idset = set([])
    for s in SeqIO.parse(filename, 'fasta'):
        oh.write('%s\t%s\n' % (s.id, str(s.seq)))
        assert s.id not in idset, 'Sequences must have all different ids'
        idset.add(s.id)
        nt_len = len(s.seq)  # save for later
    calog.info('Found %d sequences' % len(idset))
    del(idset)
    oh.close()

    if frame:
        seq_ids, seq_nts, seq_aas, stop_codons = \
            run_codon_align(table_file, frame, master_seq, fs_comp, out_dir)
    else:
        print('Looking for best frame', file=sys.stderr)
        calog.info('Running in auto frame selection')
        if nt_len < 50:
            calog.error('Sequence too short for estimating best frame')
            sys.exit('Only possible for > 50 nucleotides')

        # run alignment for all three frames
        max_len = 0
        best_frame = 123
        lens = []
        for f in [1, 2, 3]:
            seq_ids, seq_nts, seq_aas, stop_codons = \
                run_codon_align(table_file, f, master_seq, fs_comp, out_dir)
            # read first aligned aa sequence and remove '#'s
            with open('%s.aa%d' % (table_file, f), 'r') as ahin:
                first_line = ahin.readline()
            seq1 = first_line.strip().split('\t')[1].replace('#', '')
            # get maximal length of ORFs
            len1 = max([len(a) for a in seq1.split('*')])
            lens.append(len1)
            if len1 > max_len:
                max_len = len1
                best_frame = f

        # Printing information on the search for best frame
        for f in [1, 2, 3]:
            if float(lens[f - 1] - max_len) / max_len < 0.3:
                print("Frame %d has a %d bp long ORF:" %
                      (f, lens[f - 1]), end=' ', file=sys.stderr)
                if f == best_frame:
                    print("this is the best frame", file=sys.stderr)
                else:
                    print("within 30% of the best", file=sys.stderr)

        # convert the best into fasta format, nucleotides
        bn = open('codon_aligned_nt.fasta', 'w')
        with open('%s.nt%d' % (table_file, best_frame), 'r') as ahc:
            for l3 in ahc:
                cid, cseq = l3.strip().split('\t')
                bn.write('>%s\n' % cid)
                bn.write('%s\n' % fill(cseq, 80))
        bn.close()
        # and amminoacids
        ba = open('codon_aligned_aa.fasta', 'w')
        with open('%s.aa%d' % (table_file, best_frame), 'r') as ahc:
            for l3 in ahc:
                cid, cseq = l3.strip().split('\t')
                ba.write('>%s\n' % cid)
                ba.write('%s\n' % fill(cseq, 80))
        ba.close()

    # move log file and remove out_dir
    shutil.copy('%s/log.txt' % out_dir, './')
    shutil.rmtree(out_dir)


def run_codon_align(table_file, fr, ma_seq, fs, od):
    '''
    Calls the perl program codonalign
    '''
    import subprocess
    cml = '%s/codonalign %s %d %d %d %s' % (dn, table_file, ma_seq, fr, fs, od)
    ca_out = subprocess.Popen(cml, shell=True, stdout=subprocess.PIPE).stdout
    # save the output
    names = []
    nucs = []
    aas = []
    nh = open('%s.nt%d' % (table_file, fr), 'w')
    ah = open('%s.aa%d' % (table_file, fr), 'w')
    for l in ca_out:
        name, ntseq, aaseq = l.strip().split('\t')
        nh.write('%s\t%s\n' % (name, ntseq))
        ah.write('%s\t%s\n' % (name, aaseq))
        names.append(name)
        nucs.append(ntseq)
        aas.append(aaseq)
    nh.close()
    ah.close()
    # count the stop codons on the master sequence, excluding the last one
    master_stop_codons = aas[0][:-1].count('*')

    return names, nucs, aas, master_stop_codons

if __name__ == "__main__":

    import argparse
    # parse command line
    parser = argparse.ArgumentParser(description='Run codon align as in\
    http://www.hiv.lanl.gov/content/sequence/CodonAlign/codonalign.html')
    parser.add_argument('-f', '--file', dest='filename', type=str,
                        default='sample.fasta',
                        help='input file in fasta format <%(default)s>')
    parser.add_argument('-r', '--frame', dest='frame', type=int,
                        help='frame [1, 2 or 3], otherwise best is determined')
    parser.add_argument('-s', '--fs_comp', dest='fs_comp', type=int, default=5,
                        help='max n of compensating frameshifts <%(default)s>')

    args = parser.parse_args()
    main(filename=args.filename, frame=args.frame, fs_comp=args.fs_comp)
