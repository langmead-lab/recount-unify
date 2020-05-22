#!/usr/bin/env python
"""
rip_annotated_junctions.py

Rips junctions from annotation files contained in
jan_24_2016_annotations.tar.gz, as described in annotation_definition.md.
Junctions are dumped to stdout, which we record as annotated_junctions.tsv.gz
in runs/sra (same directory as this file). annotated_junctions.tsv.gz is
required by tables.py. The format of annotated_junctions.tsv.gz is

(tab-separated fields), one per junction
1. Chromosome
2. Start position (1-based, inclusive)
3. End position (1-based, inclusive)
4. Strand (+ or -)
5. anno source (abbreviation)

Must have

http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/
    hg38ToHg19.over.chain.gz

and liftOver executable available from
    https://genome-store.ucsc.edu/products/ .

Stats are written to stderr; we store them in
    annotated_junctions_stats.txt. We store hg38 regions that do not
    map to hg19 in unmapped_hg38.bed .

From the runs/sra/v2 directory, we ran

pypy rip_annotated_junctions.py
    --extract-script-dir /path/to/hisat2-2.0.1-beta
    --annotations path/to/jan_24_2016_annotations.tar.gz
    --chain /path/to/hg19ToHg38.over.chain
    --liftover /path/to/liftOver
    --unmapped unmapped_hg19.bed 2>annotated_junctions_stats.txt
    | sort -k1,1 -k2,2n -k3,3n | gzip >annotated_junctions.tsv.gz
"""

import subprocess
import tarfile
import argparse
import tempfile
import atexit
import shutil
import glob
import os
import gzip
import sys
import re

#latest gencode rel num
CURRENT_MOUSE_REL_NUM=24
CURRENT_HUMAN_REL_NUM=33

hg_file2source = {"hg19/gencode.v19.annotation.gtf.gz":"gC19","hg19/refGene.txt.gz":"rG19","hg19/acembly.txt.gz":"aC19","hg19/ccdsGene.txt.gz":"cG19","hg19/vegaGene.txt.gz":"vG19","hg19/knownGene.txt.gz":"kG19","hg19/mgcGenes.txt.gz":"mG19","hg19/lincRNAsTranscripts.txt.gz":"lR19","hg19/sibGene.txt.gz":"sG19","hg38/refGene.txt.gz":"rG38","hg38/ccdsGene.txt.gz":"cG38","hg38/knownGene.txt.gz":"kG38","hg38/mgcGenes.txt.gz":"mG38","hg38/lincRNAsTranscripts.txt.gz":"lR38","hg38/sibGene.txt.gz":"sG38","hg38/chess2.2_assembly.gtf.gz":"cH38","hg38/gencode.v24.annotation.gtf.gz":"gC24","hg38/gencode.v25.annotation.gtf.gz":"gC25","hg38/gencode.v26.annotation.gtf.gz":"gC26","hg38/gencode.v33.annotation.gtf.gz":"gC33","hg38/gencode.v29.annotation.gtf.gz":"gC29"} #,"hg38/knownAlt.txt.gz":"kA38","hg19/knownAlt.txt.gz":"kA19"}

for i in range(3,20):
    hg_file2source["hg19/gencode.v%d.annotation.gtf.gz" % (i)] = "gC%02.0f" % (i)

for i in range(20,CURRENT_HUMAN_REL_NUM+1):
    hg_file2source["hg38/gencode.v%d.annotation.gtf.gz" % (i)] = "gC%02.0f" % (i)


mg_file2source = {"m37/ccdsGene.txt.gz":"cG37","m37/refGene.txt.gz":"rG37","m37/mgcGenes.txt.gz":"mG37","m37/knownGene.txt.gz":"kG37","m37/vegaGene.txt.gz":"vG37","m37/acembly.txt.gz":"aC37","m38/GSE72311_lncrna.gtf.gz":"lR38","m38/mgcGenes.txt.gz":"mG38","m38/ccdsGene.txt.gz":"cG38","m38/knownGene.txt.gz":"kG38","m38/refGene.txt.gz":"rG38","m37/gencode.vM1.annotation.gtf.gz":"gC01","m37/fantom4_all_mrna.gtf.gz":"fM38","m37/fantom4_all_est.gtf.gz":"fE38"}
for i in range(2,CURRENT_MOUSE_REL_NUM+1):
    mg_file2source["m38/gencode.vM%d.annotation.gtf.gz" % (i)] = "gC%02.0f" % (i)

file2sources = {'hg38':hg_file2source, 'hg19':hg_file2source, 'm38':mg_file2source}

def liftover(args, annotated_junctions_from, annotated_junctions_to):
    #can handle either direction for human, but not mouse
    #if args.org == 'hg38' or args.org == 'hg19':
    if True:
        temp_from = os.path.join(extract_destination, 'from.bed')
        temp_to = os.path.join(extract_destination, 'to.bed')
        with open(temp_from, 'w') as from_stream:
            for i, junction in enumerate(annotated_junctions_from):
                # Handle incorrect junctions
                print >>from_stream, '{}\t{}\t{}\t{}\t1\t{}'.format(
                        junction[0], junction[1], junction[2], junction[4], junction[3]
                    )
        liftover_process = subprocess.check_call(' '.join([
                                                args.liftover,
                                                temp_from,
                                                args.chain,
                                                temp_to,
                                                args.unmapped
                                            ]),
                                            shell=True,
                                            executable='/bin/bash'
                                        )
        # Add all new junctions to to set
        before_liftover = len([junction for junction
                                in annotated_junctions_to
                                if junction[0] in refs])
        print >>sys.stderr, ('Below, if an RNAME is not in the chromosomal '
                             'assembly, it\'s ignored.')
        with open(temp_to) as to_stream:
            for line in to_stream:
                chrom, start, end, name, score, strand = line.strip().split('\t')[:6]
                if chrom in refs:
                    annotated_junctions_to.add((chrom, int(start), int(end), strand, name))
                else:
                    print >>sys.stderr, '({}, {}, {}, {}) not recorded.'.format(chrom, start, end, name)
        after_liftover = len([junction for junction
                                in annotated_junctions_to
                                if junction[0] in refs])
        print >>sys.stderr, ('{} annotations contributed {} junctions, and '
                             'liftover of hg19 annotations contributed an '
                             'additional {} junctions.').format(args.org,
                                    before_liftover,after_liftover - before_liftover)


if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--extract-script-dir', type=str, required=True,
            help=('path to directory containing contents of HISAT2; we '
                  'unpacked ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/'
                  'downloads/hisat-2.0.0-beta-Linux_x86_64.zip to get this')
        )
    parser.add_argument('--annotations', type=str, required=True,
            help=('annotations archive; this could be '
                  'jan_24_2016_annotations.tar.gz')
        )
    parser.add_argument('--liftover', type=str, required=False,
            help=('path to liftOver executable available from '
                  'https://genome-store.ucsc.edu/products/')
        )
    parser.add_argument('--chain', type=str, required=False,
            help=('path to unzipped liftover chain; this could be '
                  'hg38ToHg19.over.chain')
        )
    parser.add_argument('--unmapped', type=str, required=False,
            help='BED in which unmapped junctions should be stored'
        )
    parser.add_argument('--org', type=str, required=True,
            help='Organism reference version (hg38, hg19, m38)'
        )
    parser.add_argument('--org-sizes', type=str, required=True,
            help='Organism reference version (hg38, hg19, m38) chromosome sizes file'
        )
    parser.add_argument('--temp', type=str, default='./tmp',
            help='temporary directory, [./tmp]'
        )
    args = parser.parse_args()
    file2source = file2sources[args.org]    
    #extract_destination = tempfile.mkdtemp()
    extract_destination = args.temp
    #atexit.register(shutil.rmtree, extract_destination)
    with tarfile.open(args.annotations, 'r:gz') as tar:
        tar.extractall(path=extract_destination)
    extract_splice_sites_path = os.path.join(args.extract_script_dir,
                                                'extract_splice_sites.py')
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    with open(args.org_sizes) as hg38_stream:
        refs = set([tokens.strip().split('\t')[0] for tokens in hg38_stream])
    annotated_junctions_hg19 = set()
    annotated_junctions_hg38 = set()
    ref_patt = re.compile(r'\/(((hg38)|(m38)|(hg19)|(m37))\/)')
    for junction_file in glob.glob(
            os.path.join(extract_destination, 'anno', '*', '*')):
        sys.stderr.write("about to extract jx from %s\n" % (junction_file))
        m = ref_patt.search(junction_file)
        if m is None:
            sys.stderr.write("junction file with unknown reference label %s, terminating\n" % junction_file)
            sys.exit(-1)
        label = m.group(1) + os.path.basename(junction_file)
        datasource_code = file2source[label]
        unique_junctions = set()
        if 'hg38' in junction_file or 'm38' in junction_file:
            annotated_junctions = annotated_junctions_hg38
        else:
            annotated_junctions = annotated_junctions_hg19
        if 'gencode' in junction_file or 'chess' in junction_file or 'GSE72311' in junction_file or 'fantom' in junction_file:
            #extract_splice_sites_path prints 0-based, exon coords around junctions
            #hence the +2 for the start here
            extract_process = subprocess.Popen(' '.join([sys.executable,
                                            extract_splice_sites_path,
                                            '<(gzip -cd %s) > %s.extracted'
                                               % (junction_file,junction_file)]),
                                        shell=True,
                                        executable='/bin/bash',
                                        stdout=subprocess.PIPE)
            exit_code = extract_process.wait()
            if exit_code != 0:
                raise RuntimeError(
                    'extract_splice_sites.py had nonzero exit code {}.'.format(exit_code))
            fin = open("%s.extracted" % (junction_file), "r")
            for line in fin:
                tokens = line.strip().split('\t')
                tokens[1] = int(tokens[1]) + 2
                tokens[2] = int(tokens[2])
                if tokens[2] < tokens[1]:
                    print >>sys.stderr, (
                            'Invalid junction ({}, {}, {}) from file {}. '
                            'Skipping.').format(tokens[0], tokens[1], tokens[2], junction_file)
                    continue
                tokens[-1] = datasource_code + ':::' + tokens[-1]
                junction_to_add = tuple(tokens)
                annotated_junctions.add(junction_to_add)
                unique_junctions.add(junction_to_add)
            fin.close()
        else:
            transcript_id_col = 1
            offset = 1
            if 'knownGene' in junction_file:
                offset = 0
                transcript_id_col = 11
            for line in gzip.open(junction_file):
                line = line.rstrip()
                tokens = line.split('\t')
                exons = [(int(start), int(end)) for start, end in
                                zip(tokens[8+offset].split(','),
                                    tokens[9+offset].split(','))[:-1]]
                junctions_to_add = [(tokens[1+offset], exons[i-1][1] + 1, exons[i][0],
                         tokens[2+offset], "%s:::%s" % (datasource_code, tokens[transcript_id_col]))
                        for i in xrange(1, len(exons))]
                final_junctions_to_add = []
                for junction in junctions_to_add:
                    if junction[2] < junction[1]:
                        print >>sys.stderr, (
                                'Invalid junction ({}, {}, {}) from line {} from file {}. '
                                'Skipping.').format(
                                    junction[0], junction[1], junction[2], line,
                                    junction_file)
                        continue
                    final_junctions_to_add.append(junction)
                annotated_junctions.update(final_junctions_to_add)
                unique_junctions.update(final_junctions_to_add)
        print >>sys.stderr, 'Junctions in {}: {}'.format(
                label,
                len(unique_junctions)
            )
    #now do liftover---either direction depending on what's passed in via args.org
    annotated_junctions = annotated_junctions_hg38
    if args.org == 'hg38' or args.org == 'm38':
        liftover(args, annotated_junctions_hg19, annotated_junctions_hg38)
    elif args.org == 'hg19':
        liftover(args, annotated_junctions_hg38, annotated_junctions_hg19)
        annotated_junctions = annotated_junctions_hg19
    junc2datasource = {}
    junc2transcript = {}
    for junction in annotated_junctions:
        if junction[0] in refs:
            (source, transcript) = re.split(':::',junction[4])
            if junction[:4] not in junc2datasource:
                junc2datasource[junction[:4]]=set()
            junc2datasource[junction[:4]].add(source)
            if junction[:4] not in junc2transcript:
                junc2transcript[junction[:4]]=set()
            junc2transcript[junction[:4]].add(junction[4])
    seen = set()
    for junction in annotated_junctions:
        if junction[0] in refs and junction[:4] not in seen:
            sources = ",".join(sorted(junc2datasource[junction[:4]]))
            transcripts = ",".join(sorted(junc2transcript[junction[:4]]))
            print "%s\t%s\t%s" % ('\t'.join(map(str, junction[:4])), sources, transcripts)
            seen.add(junction[:4])
