#!/usr/bin/env python3

import gzip
import json
import os
import re
import sys
import argparse

def create_annotation_lookup(variant_catalogue):
    """Create a lookup dictionary from the variant catalog"""
    with open(variant_catalogue, 'r') as f:
        catalog = json.load(f)

    # Create lookup dictionary by LocusId
    lookup = {}
    for entry in catalog:
        locus_id = entry['LocusId']
        lookup[locus_id] = {
            'DisplayRU': entry.get('DisplayRU', '.'),
            'HGNCId': str(entry.get('HGNCId', '.')),
            'Disease': entry.get('Disease', '.'),
            'InheritanceMode': entry.get('InheritanceMode', '.'),
            'STR_NORMAL_MAX': str(entry.get('NormalMax', '.')),
            'STR_PATHOLOGIC_MIN': str(entry.get('PathologicMin', '.'))
        }

    print(f"Created lookup for {len(lookup)} loci")
    return lookup

def annotate_vcf(vcf_file, lookup_data, sample_id):
    """Annotate VCF file with STR information"""
    annotated_count = 0
    total_variants = 0
    header_processed = False

    with gzip.open(vcf_file, 'rt') as infile, open(f'{sample_id}.vcf', 'w') as outfile:
        for line in infile:
            if line.startswith('##fileformat='):
                # Change VCF version to 4.2
                outfile.write('##fileformat=VCFv4.2\n')
            elif line.startswith('##INFO=<ID=LOCUS,'):
                # Skip the original LOCUS line - we'll replace it with LocusId
                continue
            elif line.startswith('##INFO=<ID=LocusId,'):
                # Write the LocusId line and add new annotation fields
                outfile.write('##INFO=<ID=LocusId,Number=1,Type=String,Description="Locus ID">\n')
                outfile.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit (motif) for STRs">\n')
                outfile.write('##INFO=<ID=HGNCId,Number=1,Type=Integer,Description="HGNC gene id for associated disease gene">\n')
                outfile.write('##INFO=<ID=Disease,Number=1,Type=String,Description="Associated disorder">\n')
                outfile.write('##INFO=<ID=InheritanceMode,Number=1,Type=String,Description="Main mode of inheritance for disorder">\n')
                outfile.write('##INFO=<ID=STR_NORMAL_MAX,Number=1,Type=Integer,Description="Max number of repeats allowed to call as normal">\n')
                outfile.write('##INFO=<ID=STR_PATHOLOGIC_MIN,Number=1,Type=Integer,Description="Min number of repeats required to call as pathologic">\n')
                header_processed = True
            elif line.startswith('##') and not header_processed and 'INFO=<ID=RUS_REF' in line:
                # If we haven't processed the header yet and we're at RUS_REF, add our fields before it
                outfile.write('##INFO=<ID=LocusId,Number=1,Type=String,Description="Locus ID">\n')
                outfile.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit (motif) for STRs">\n')
                outfile.write('##INFO=<ID=HGNCId,Number=1,Type=Integer,Description="HGNC gene id for associated disease gene">\n')
                outfile.write('##INFO=<ID=Disease,Number=1,Type=String,Description="Associated disorder">\n')
                outfile.write('##INFO=<ID=InheritanceMode,Number=1,Type=String,Description="Main mode of inheritance for disorder">\n')
                outfile.write('##INFO=<ID=STR_NORMAL_MAX,Number=1,Type=Integer,Description="Max number of repeats allowed to call as normal">\n')
                outfile.write('##INFO=<ID=STR_PATHOLOGIC_MIN,Number=1,Type=Integer,Description="Min number of repeats required to call as pathologic">\n')
                outfile.write(line)
                header_processed = True
            elif line.startswith('##'):
                # Write other header lines as-is
                outfile.write(line)
            elif line.startswith('#CHROM'):
                # Write column header
                outfile.write(line)
            else:
                # Process variant lines
                total_variants += 1
                fields = line.strip().split('\t')
                info_field = fields[7]

                # Replace LOCUS= with LocusId= in INFO field
                info_field = info_field.replace('LOCUS=', 'LocusId=')

                # Extract LocusId value
                locus_match = re.search(r'LocusId=([^;]+)', info_field)
                if locus_match:
                    locus_id = locus_match.group(1)
                    fields[2] = locus_id  # Set ID field to LocusId value

                    # Add annotations if available
                    if locus_id in lookup_data:
                        ann_data = lookup_data[locus_id]
                        new_annotations = []

                        for key, value in ann_data.items():
                            if value != '.':
                                if key == 'DisplayRU':
                                    new_annotations.append(f"RU={value}")
                                else:
                                    new_annotations.append(f"{key}={value}")

                        if new_annotations:
                            info_field += ';' + ';'.join(new_annotations)
                            annotated_count += 1

                fields[7] = info_field
                outfile.write('\t'.join(fields) + '\n')
    print(f"Writing annotated VCF to: {os.path.abspath(f'{sample_id}.vcf')}")
    print(f"Processed {total_variants} variants, annotated {annotated_count}")
    return total_variants, annotated_count

def create_tsv_summary(sample_id):
    """Create TSV summary from annotated VCF"""
    with open(f'{sample_id}_annotated.tsv', 'w') as tsv_file:
        # Write header
        tsv_file.write("CHROM\tPOS\tID\tREF\tALT\tLocusId\tRU\tHGNCId\tDisease\tInheritanceMode\tSTR_NORMAL_MAX\tSTR_PATHOLOGIC_MIN\n")

        # Process VCF lines
        with gzip.open(f'{sample_id}.vcf.gz', 'rt') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                chrom, pos, id_field, ref, alt = fields[:5]
                info = fields[7]

                # Extract fields from INFO
                info_dict = {}
                for part in info.split(';'):
                    if '=' in part:
                        key, value = part.split('=', 1)
                        info_dict[key] = value

                # Get annotation values
                locus_id = info_dict.get('LocusId', '')
                ru = info_dict.get('RU', '')
                hgnc_id = info_dict.get('HGNCId', '')
                disease = info_dict.get('Disease', '')
                inheritance = info_dict.get('InheritanceMode', '')
                normal_max = info_dict.get('STR_NORMAL_MAX', '')
                pathologic_min = info_dict.get('STR_PATHOLOGIC_MIN', '')

                tsv_file.write(f"{chrom}\t{pos}\t{id_field}\t{ref}\t{alt}\t{locus_id}\t{ru}\t{hgnc_id}\t{disease}\t{inheritance}\t{normal_max}\t{pathologic_min}\n")

def main():
    parser = argparse.ArgumentParser(description='Annotate STR VCF files with disease information')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--catalogue', required=True, help='Variant catalogue JSON file')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    parser.add_argument('--output-dir', default='.', help='Output directory')

    args = parser.parse_args()

    # Create annotation lookup
    lookup_data = create_annotation_lookup(args.catalogue)

    # Annotate VCF
    total_variants, annotated_count = annotate_vcf(args.vcf, lookup_data, args.sample_id)

    # Write log file
    with open(f'{args.sample_id}_annotation.log', 'w') as log:
        log.write(f"Sample: {args.sample_id}\n")
        log.write(f"Total variants: {total_variants}\n")
        log.write(f"Annotated variants: {annotated_count}\n")
        if total_variants > 0:
            log.write(f"Annotation rate: {annotated_count/total_variants*100:.1f}%\n")

if __name__ == '__main__':
    main()
