"""Identify protein coding genes from GTF reference file for GTEx."""

import argparse
import os
from collections import namedtuple

parser = argparse.ArgumentParser(prog=__file__, description=__doc__)
parser.add_argument('gtf_path')
parser.add_argument('-o', '--output_path', default='./')

GTF_fields = ['seqname', 'source', 'feature',
              'start', 'end',
              'score',
              'strand', 'frame',
              'attribute']
GTFRecord = namedtuple('GTFRecord', GTF_fields)

output_fields = ['gene_id', 'gene_type', 'gene_name']
OutputRecord = namedtuple('OutputRecord', output_fields)

if __name__ == '__main__':
    args = parser.parse_args()

    records = []
    with open(args.gtf_path) as file:
        for line in file:
            if line.startswith('##'):
                continue

            gtf_record = GTFRecord(*line.rstrip(';\n').split('\t'))
            if gtf_record.feature != 'gene':  # Skip all non-gene entries
                continue

            attributes = {}
            for attribute in gtf_record.attribute.split(';'):
                tag, value = attribute.lstrip().rstrip().split(' ')
                attributes[tag] = value.lstrip('"').rstrip('"')

            records.append(OutputRecord(*[attributes[field_name] for field_name in output_fields]))

    prefix = args.output_path
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    with open(f'{prefix}/genes.tsv', 'w') as file:
        file.write('\t'.join(output_fields) + '\n')
        for record in records:
            file.write('\t'.join(record) + '\n')
