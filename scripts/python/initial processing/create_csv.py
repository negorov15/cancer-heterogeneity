# Create a CSV file with cell line information, including barcodes and mutated genes, for downstream processing. 
# The output file is saved as 'resources/cell_lines.tsv'.

import csv
from pathlib import Path

header = ['cell_line', 'pool', 'barcodes', 'mutated_genes']
rows = [
    ['MDA-MB-231', 'pool2', 'resources/barcodes/mdamb231_barcodes.txt', 'resources/gene_list_bed/MDA-MB-231_without_expressions_mutated.bed'],
    ['SCC-4',      'pool2', 'resources/barcodes/scc4_barcodes.txt',     'resources/gene_list_bed/SCC-4_without_expressions_mutated'],
    ['A-253',      'pool3', 'resources/barcodes/a253_barcodes.txt',     'resources/gene_list_bed/A-253_without_expressions_mutated.bed'],
    ['MCF-7',      'pool3', 'resources/barcodes/7mcf_barcodes.txt',     'resources/gene_list_bed/MCF7_without_expressions_mutated'],
    ['HCT-116',    'pool3', 'resources/barcodes/hct116_barcodes.txt',   'resources/gene_list_bed/HCT116_without_expressions_mutated.bed'],
    ['BT-474',     'pool4', 'resources/barcodes/bt474_barcodes.txt',    'resources/gene_list_bed/BT-474_without_expressions_mutated.bed'],
    ['HCT-15',     'pool4', 'resources/barcodes/hct15_barcodes.txt',    'resources/gene_list_bed/HCT-15_without_expressions_mutated.bed'],
]

out = Path('resources/cell_lines.tsv')
out.parent.mkdir(parents=True, exist_ok=True)

with out.open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(header)
    w.writerows(rows)

# data = [
#     ['Cell line', 'Model_ID', 'path_to_count_matrix'],
#     ['Caco-2', 'ACH-000003', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/Caco-2.count_mtx.tsv.txt'],
#     ['MDA-MB-231', 'ACH-000768', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/MDA-MB-231.count_mtx.tsv.txt'],
#     ['SCC-4', 'ACH-000238', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/SCC-4.count_mtx.tsv.txt'],
#     ['A-253', 'ACH-000740', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/A-253.count_mtx.tsv.txt'],
#     ['MCF-7', 'ACH-000019', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/MCF7.count_mtx.tsv.txt'],
#     ['HCT-116', 'ACH-000971', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/HCT116.count_mtx.tsv.txt'],
#     ['BT-474', 'ACH-000927', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/BT-474.count_mtx.tsv.txt'],
#     ['HCT-15', 'ACH-000997', '/Users/nikoegorov/Documents/internship/data/CSE0000121_processed/HCT-15.count_mtx.tsv.txt'],
# ]
