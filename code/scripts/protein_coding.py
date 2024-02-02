"""Identify protein coding genes from GTF reference file for GTEx."""

import os
from collections import namedtuple

GTF_fields = ['seqname', 'source', 'feature',
              'start', 'end',
              'score',
              'strand', 'frame',
              'attribute']
GTFRecord = namedtuple('GTFRecord', GTF_fields)

output_fields = ['gene_id', 'gene_type', 'gene_name']
OutputRecord = namedtuple('OutputRecord', output_fields)

data_path = '../../data/GTEx/gene_model/references_v8_gencode.v26.GRCh38.genes.gtf'

records = []
with open(data_path) as file:
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

if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/genes.tsv', 'w') as file:
    file.write('\t'.join(output_fields) + '\n')
    for record in records:
        file.write('\t'.join(record) + '\n')
