from Bio import SeqIO
import csv
import os
import gzip
import re
from mimetypes import guess_type
from functools import partial
import argparse
import itertools
import sys
import pathlib

barcode_re = re.compile(":[acgtn]+$", re.IGNORECASE)

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("barcode_map", type=argparse.FileType('r'),help="Tab delimited file giving barcodes and the corresponding sample")
  parser.add_argument("fasta_output", type=argparse.FileType('a'),help="Path for FASTA output")
  parser.add_argument("--fastqdir", help="Directory containing FASTQ files", default=None)
  parser.add_argument('--debug', type=int, default=None, help='Number of reads to parse from each FASTQ')  
  args = parser.parse_args()
  
  if args.fastqdir:
    fastq_dir = args.fastqdir
  else:
    fastq_dir = os.path.dirname(args.barcode_map.name)


    
      

  barcode_dict = parseBarcodeMap(args.barcode_map)
  # for x in barcode_dict:
  #   print(x)
  args.fasta_output.truncate(0) # to avoid adding

  

  
  #----------------------------------------------------------------------------
  # out_path = pathlib.PurePath(args.fasta_output.name).parent
  # print(out_path)
  
  generated_mappath = pathlib.Path(args.barcode_map.name).stem + "_to2column.csv"
  generated_mappath = pathlib.Path(args.fasta_output.name).with_name(generated_mappath)
  print("New Map File:", generated_mappath)
  # sys.exit(1)
  with generated_mappath.open('w') as csvfile:
    mapwriter = csv.writer(csvfile, delimiter='\t',
                             quoting=csv.QUOTE_MINIMAL)
    for barcode, samplename, fastq_filename in barcode_dict:
      print(barcode, samplename, fastq_filename)
      mapwriter.writerow((barcode, samplename))
      fastq_path = os.path.join(fastq_dir, fastq_filename)
      fastqToFASTA(fastq_path, barcode, args.fasta_output, args.debug)
  #----------------------------------------------------------------------------
  
  

def parseBarcodeMap(csvfile):
  barcode_map_reader = csv.reader(csvfile, delimiter='\t')
  for row in barcode_map_reader:
    yield(row)

def fastqToFASTA(fastq_filename, corrected_barcode=None, fasta_handle=None, debug=None):
  readcount=0
  encoding = guess_type(fastq_filename)[1]  # uses file extension
  if encoding is None:
    _open = open
  elif encoding == 'gzip':
    _open = partial(gzip.open, mode='rt')
  else:
    raise ValueError('Unknown file encoding: "{}"'.format(encoding))
  with _open(fastq_filename) as f:
    for record in itertools.islice(SeqIO.parse(f, 'fastq'), debug):
        barcode_replaced = barcode_re.sub(":" + corrected_barcode, record.description)
        record.description = barcode_replaced
        SeqIO.write(record, fasta_handle, "fasta")
    if (debug):
      print("DEBUG MODE: processed {0} reads".format(debug), file=sys.stderr)

if __name__ == '__main__':
  main()
