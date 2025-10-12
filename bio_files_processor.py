"""
bio_files_processor.py

Utilities for processing BLAST, FASTA, and GenBank files:
- convert_multiline_fasta_to_oneline: converts multi-line FASTA to one-line format
- parse_blast_output: extracts top hits from BLAST text output
- select_genes_from_gbk_to_fasta: extracts neighboring genes from GenBank files
"""

import os
from typing import Union, List, Set


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str = None
) -> str:
    """
    Converts a multi-line FASTA file to one-line format.
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output file (optional). If None, uses input name with '_oneline' suffix
    
    """
    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = f"{base}_oneline{ext if ext else '.fasta'}"
    
    with open(input_fasta, 'r', encoding='utf-8') as fin, \
         open(output_fasta, 'w', encoding='utf-8') as fout:
        
        current_header = None
        current_seq = []
        
        for line in fin:
            line = line.rstrip('\r\n')
            
            if line.startswith('>'):
                if current_header is not None:
                    fout.write(current_header + '\n')
                    fout.write(''.join(current_seq) + '\n')
                
                current_header = line
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # Write last sequence
        if current_header is not None:
            fout.write(current_header + '\n')
            fout.write(''.join(current_seq) + '\n')
    
    return output_fasta


def parse_blast_output(input_file: str, output_file: str) -> str:
    """
    Parses BLAST text output and extracts the top hit description for each query.
    
    Args:
        input_file: Path to BLAST output text file
        output_file: Path to output file with top descriptions
    """
    descriptions = []
    in_alignments_section = False
    found_first_hit = False
    desc_col_end = None

    with open(input_file, 'r', encoding='utf-8') as fin:
        for line in fin:
            line = line.rstrip('\r\n')

            if 'Sequences producing significant alignments:' in line:
                in_alignments_section = True
                found_first_hit = False
                desc_col_end = None
                continue

            if in_alignments_section:
                # Find header line with "Description" to determine column width
                if 'Description' in line and desc_col_end is None:
                    for pattern in ['Max Score', 'Total Score', 'Query Cover', 'E value', 'Per.']:
                        idx = line.find(pattern)
                        if idx != -1:
                            desc_col_end = idx
                            break

                    if desc_col_end is None:
                        desc_start = line.find('Description')
                        if desc_start != -1:
                            search_start = desc_start + len('Description')
                            for i in range(search_start, len(line) - 1):
                                if line[i:i+2] == '  ' and line[i+2:i+4].strip():
                                    desc_col_end = i
                                    break
                    continue

                if line.startswith('---'):
                    continue

                if not line.strip():
                    if found_first_hit:
                        in_alignments_section = False
                    continue

                if not found_first_hit and desc_col_end is not None:
                    description = line[:desc_col_end].strip()
                    if description:
                        descriptions.append(description)
                        found_first_hit = True
                        in_alignments_section = False

    if not descriptions:
        raise ValueError("No hit descriptions found in BLAST output.")

    with open(output_file, 'w', encoding='utf-8') as fout:
        for desc in descriptions:
            fout.write(desc + '\n')

    return output_file



def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: Union[str, List[str]],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = 'selected_genes.fasta'
) -> str:
    """
    Extracts neighboring genes from a GenBank file based on genes of interest.
    
    Args:
        input_gbk: Path to input GenBank file
        genes: Gene name(s) of interest (string or list of strings)
        n_before: Number of genes before each target gene (default: 1)
        n_after: Number of genes after each target gene (default: 1)
        output_fasta: Path to output FASTA file
    
    """
    if isinstance(genes, str):
        genes = [genes]
    genes_set = set(genes)
    
    if n_before < 1 or n_after < 1:
        raise ValueError("n_before and n_after must be >= 1")
    
    # Parse GenBank file and extract all CDS features with translations
    all_genes = []
    
    with open(input_gbk, 'r', encoding='utf-8') as fin:
        in_cds = False
        current_gene = {}
        current_translation = []
        
        for line in fin:
            line_stripped = line.rstrip('\r\n')
            num_spaces = len(line) - len(line.lstrip(' '))
            
            # Detect start of CDS feature (5 spaces + "CDS")
            if num_spaces == 5 and line_stripped.strip().startswith('CDS'):
                # Save previous CDS if it has translation
                if current_translation:
                    trans = ''.join(current_translation)
                    all_genes.append({
                        'gene': current_gene.get('gene'),
                        'locus_tag': current_gene.get('locus_tag'),
                        'product': current_gene.get('product'),
                        'translation': trans
                    })
                
                in_cds = True
                current_gene = {'gene': None, 'locus_tag': None, 'product': None}
                current_translation = []
                continue
            
            if in_cds:
                # Check for new feature 
                if num_spaces <= 5 and line_stripped.strip():
                    if current_translation:
                        trans = ''.join(current_translation)
                        all_genes.append({
                            'gene': current_gene.get('gene'),
                            'locus_tag': current_gene.get('locus_tag'),
                            'product': current_gene.get('product'),
                            'translation': trans
                        })
                    in_cds = False
                    current_gene = {}
                    current_translation = []
                    continue
                
                if '/gene=' in line_stripped:
                    gene_name = line_stripped.split('/gene=')[1].strip('"')
                    current_gene['gene'] = gene_name
                
                if '/locus_tag=' in line_stripped:
                    locus = line_stripped.split('/locus_tag=')[1].strip('"')
                    current_gene['locus_tag'] = locus
                
                if '/product=' in line_stripped:
                    product = line_stripped.split('/product=')[1].strip('"')
                    current_gene['product'] = product
                
                if '/translation=' in line_stripped:
                    trans_start = line_stripped.split('/translation=')[1].strip('"')
                    current_translation.append(trans_start)
                elif current_translation and num_spaces == 21:
                    # Continue multi-line translation (only lines with 21 spaces)
                    current_translation.append(line_stripped.strip().strip('"'))
        
        # Save last CDS if exists
        if current_translation:
            trans = ''.join(current_translation)
            all_genes.append({
                'gene': current_gene.get('gene'),
                'locus_tag': current_gene.get('locus_tag'),
                'product': current_gene.get('product'),
                'translation': trans
            })

    target_indices = []
    for i, gene_info in enumerate(all_genes):
        gene_name = gene_info.get('gene')
        locus_tag = gene_info.get('locus_tag')
        
        if gene_name in genes_set or locus_tag in genes_set:
            target_indices.append(i)
    
    # Collect neighboring genes 
    neighbors_indices: Set[int] = set()
    
    for idx in target_indices:
        for i in range(max(0, idx - n_before), idx):
            neighbors_indices.add(i)
        
        for i in range(idx + 1, min(len(all_genes), idx + 1 + n_after)):
            neighbors_indices.add(i)
    
    with open(output_fasta, 'w', encoding='utf-8') as fout:
        for idx in sorted(neighbors_indices):
            gene_info = all_genes[idx]
            gene_name = gene_info.get('gene') or gene_info.get('locus_tag') or f'gene_{idx}'
            product = gene_info.get('product') or 'unknown product'
            translation = gene_info.get('translation', '')
            
            if not translation:
                continue
            
            header = f">{gene_name} {product}"
            fout.write(header + '\n')
            fout.write(translation + '\n')
    
    return output_fasta