from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import time
import logging

__author__ = 'Martin Aryee'

logger = logging.getLogger('root')

#args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
#base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
#args['read1'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
#args['read2'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
#args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
#args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')

def fq(file, start=0, max_count=-1):
    count=0
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        if start>0:
            junk = f.readlines(4*start)

        while True:
            if max_count >0 and count >max_count: break
            count+=1

            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]


def get_sample_id(i1, i2, sample_names):
    seq1 = i1[1]
    seq2 = i2[1]
    sample_barcode = seq1[1:8] + seq2[1:8]
    if sample_barcode in sample_names:
        return sample_names[sample_barcode]
    else:
        return sample_barcode


def read_core(read1,read2,index1,index2,sample_names,start_record,stride):
        r1s = [fq(read1,start=start_record,count=stride)]
        r2s = [fq(read2,start=start_record,count=stride)]
        i1s = [fq(index1,start=start_record,count=stride)]
        i2s = [fq(index2,start=start_record,count=stride)]

        if len(r1s) == 0:
            return None

        ids = [get_sample_id(i1,i2,sample_names) for i1, i2 in i1s[start_record:end_record], i2s[start_record:end_record]]

        keys = set(ids)
        r1_map = dict([(k,[]) for k in keys)
        r2_map = dict([(k,[]) for k in keys)
        i1_map = dict([(k,[]) for k in keys)
        i2_map = dict([(k,[]) for k in keys)

        for i,e in ids:
            r1_map[ids[i]] = r1s[i]
            r2_map[ids[i]] = r2s[i]
            i1_map[ids[i]] = i1s[i]
            i2_map[ids[i]] = i2s[i]

        return r1_map, r2_map, i1_map, i2_map


def demultiplex(read1, read2, index1, index2, sample_barcodes, out_dir, min_reads=10000):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if type(sample_barcodes) != dict:
        sample_names = {}
        if not sample_barcodes==None:
            for line in open(sample_barcodes, 'r'):
                fields = line.strip().split('\t')
                if len(fields)==2:
                    sampleid, barcode = fields
                    sample_names[barcode] = sampleid
    else:
        sample_names = sample_barcodes

    outfiles_r1 = {}
    outfiles_r2 = {}
    outfiles_i1 = {}
    outfiles_i2 = {}

    total_count = 0
    count = {}
    buffer_r1 = {}
    buffer_r2 = {}
    buffer_i1 = {}
    buffer_i2 = {}

    #it = itertools.izip(fq(args['read1']), fq(args['read2']), fq(args['index1']), fq(args['index2']))
    #for r1,r2,i1,i2 in itertools.islice(it, 0, 100):
    start = time.time()

    cores = 10
    stride = 10000
    total_count = 0

    all_r1s = {}
    all_r2s = {}
    all_i1s = {}
    all_i2s = {}

    while True
        for i,c in range(cores):

            start_record = i*stride
            end_record=(i+1)*stride
            out = read_core(read1,read2,index1,index2,sample_names,start_record,stride)
            if out == None:
                break
            else:
                total_count += len(out[0])

            all_r1s.extend(out[0])
            all_r2s.extend(out[1])
            all_i1s.extend(out[2])
            all_i2s.extend(out[3])
        break

    for sample_id in all_r1s.keys():

        if len(all_r1s[k]) >= min_reads:
            outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % sample_id), 'w')
            outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % sample_id), 'w')
            outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % sample_id), 'w')
            outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % sample_id), 'w')
            # Spill the buffers to sample-specific fastqs
            for record in buffer_r1[sample_id] + r1:
                outfiles_r1[sample_id].write(''.join(record))
            for record in buffer_r2[sample_id] + r2:
                outfiles_r2[sample_id].write(''.join(record))
            for record in buffer_i1[sample_id] + i1:
                outfiles_i1[sample_id].write(''.join(record))
            for record in buffer_i2[sample_id] + i2:
                outfiles_i2[sample_id].write(''.join(record))

            del all_r1s[sample_id]
            del all_r2s[sample_id]
            del all_i1s[sample_id]
            del all_i2s[sample_id]

    # Write remaining buffered reads to a single fastq.
    # (These reads correspond to barcodes that were seen less than min_reads times)
    undetermined_r1 = open(os.path.join(out_dir, 'undetermined.r1.fastq'), 'w')
    undetermined_r2 = open(os.path.join(out_dir, 'undetermined.r2.fastq'), 'w')
    undetermined_i1 = open(os.path.join(out_dir, 'undetermined.i1.fastq'), 'w')
    undetermined_i2 = open(os.path.join(out_dir, 'undetermined.i2.fastq'), 'w')
    for sample_id in all_r1s.keys():
        for record in all_r1s[sample_id]:
            undetermined_r1.write(''.join(record))
        for record in all_r2s[sample_id]:
            undetermined_r2.write(''.join(record))
        for record in all_i1s[sample_id]:
            undetermined_i1.write(''.join(record))
        for record in all_i2s[sample_id]:
            undetermined_i2.write(''.join(record))

    # Close files
    for sample_id in outfiles_r1:
        outfiles_r1[sample_id].close()
        outfiles_r2[sample_id].close()
        outfiles_i1[sample_id].close()
        outfiles_i2[sample_id].close()
    undetermined_r1.close()
    undetermined_r2.close()
    undetermined_i1.close()
    undetermined_i2.close()


    logger.info('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads.', len(outfiles_r1), total_count, min_reads)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1', required=True)
    parser.add_argument('--read2', required=True)
    parser.add_argument('--index1', required=True)
    parser.add_argument('--index2', required=True)
    parser.add_argument('--min_reads', type=int, default=10000)
    parser.add_argument('--sample_barcodes')
    parser.add_argument('--out_dir', default='.')
    args = vars(parser.parse_args())

    demultiplex(args['read1'], args['read2'], args['index1'], args['index2'], args['sample_barcodes'], args['out_dir'], min_reads=args['min_reads'])

if __name__ == '__main__':
    main()
