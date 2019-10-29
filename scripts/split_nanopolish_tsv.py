#!/usr/bin/env python3
# Function to split large tsv by (consecutive) contigs
# Usage: python split_nanopolish_tsv.py <nanopolish.tsv> [chunk_size:default=1e7]

import csv
import os
import sys


def split_tsv(file_path, chunk_size=1e7):
    file_name = os.path.splitext(file_path)[0]

    with open(file_path, 'r', newline='', encoding='utf-8') as fp:
        outfile = None
        writer = None
        contig = None
        parts = dict()
        lines = dict()

        reader = csv.reader(fp, delimiter='\t', quotechar='"')
        header = next(reader)

        for row in reader:
            if contig is None:
                # first line
                contig = row[0]
                parts[contig] = 0
                lines[contig] = 0

                # open first line
                outname = '{0}.{1}.{2}.tsv'.format(file_name, contig, parts[contig])
                outfile = open(outname, 'w', newline='', encoding='utf-8')
                writer = csv.writer(outfile, delimiter='\t', quotechar='"')
                writer.writerow(header)

            new_contig = row[0]

            if new_contig != contig:
                # close and start new file
                outfile.close()

                if not new_contig in parts:
                    parts[new_contig] = 0
                    lines[new_contig] = 0

                outname = '{0}.{1}.{2}.tsv'.format(file_name, new_contig, parts[new_contig])
                outfile = open(outname, 'a+', newline='', encoding='utf-8')
                writer = csv.writer(outfile, delimiter='\t', quotechar='"')
                if outfile.tell() == 0:
                    # writing new
                    writer.writerow(header)

                contig = new_contig

            elif lines[contig] >= chunk_size:
                # close and start new part
                outfile.close()
                print('File {} completes'.format(outname))

                parts[contig] += 1
                lines[contig] = 0

                outname = '{0}.{1}.{2}.tsv'.format(file_name, contig, parts[contig])
                outfile = open(outname, 'w', newline='', encoding='utf-8')
                writer = csv.writer(outfile, delimiter='\t', quotechar='"')
                writer.writerow(header)

            writer.writerow(row)
            lines[contig] += 1

        outfile.close()

if __name__ == '__main__':
    try:
        file_path = sys.argv[1]
        if len(sys.argv) == 3:
            chunk_size = int(sys.argv[2])
            split_tsv(file_path, chunk_size)
        else:
            split_tsv(file_path)
    except:
        raise Exception("Please supply correct argument(s)")
