#!/usr/bin/env python3
#python/3.7.0
import os
import sys
import glob
import gzip
import string

path2in = '/u/scratch/t/tianhao/NovaSeq100622/data/'
path2out = '/u/scratch/t/tianhao/NovaSeq100622/split/'

def split(infile1):
  outfile = infile1.rsplit('/')[-1].rsplit('_001.fastq')[0]
  #lane = infile1.rsplit('/')[-3][-1]
  chunksize = 5000000
  fid = 1
  if os.path.exists(path2out+outfile+'_001.fastq'): 
    fid = 900
  with gzip.open(infile1) as infile:
    f = open(path2out+outfile+'_%03d.fastq' %fid, 'wb')
    for i,line in enumerate(infile):
      f.write(line)
      if not (i+1)%chunksize:
        print("file%d" %fid)
        f.close()
        fid += 1
        #if fid > 500: break
        f = open(path2out+outfile+'_%03d.fastq' %fid, 'wb')
    f.close()

def main():
  infiles = sorted(glob.glob(path2in+'*/*.fastq.gz'))
  sampledict = {}
  for infile in infiles:
    print(infile)
    split(infile)

if __name__ == '__main__':
  main()


