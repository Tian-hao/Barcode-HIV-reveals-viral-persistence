#!/usr/bin/env python3

import glob
import time
from collections import Counter

workpath  = '/u/scratch/t/tianhao/NovaSeq040122/'

def main():
  infiles  = sorted(glob.glob(workpath+'ISS_RNA/RNA_barcode_*.txt'))
  samplelist = [x.rsplit('_')[-2] for x in infiles]
  samplelist = sorted(list(set(samplelist)))
  for sample in samplelist:
    summarize(sample)

def summarize(sample):
  infiles   = sorted(glob.glob(workpath+'ISS_RNA/RNA_barcode_'+sample+'_*.txt'))
  countdict = {} #countdict[UMI] = [bc1,bc2,...]
  bcdepthdict = {} #bcdepthdict[bc] = N
  start_time = time.time()
  for infile in infiles:
    inhandle = open(infile)
    print('processing '+infile+'\ntime: '+str(time.time()-start_time))
    for line in inhandle:
      line = line.rstrip().rsplit('\t')
      if len(line) < 2: continue
      UMI  = line[0]
      bc   = line[1]
      if UMI not in countdict:
        countdict[UMI] = []
      countdict[UMI].append(bc)
      if bc not in bcdepthdict:
        bcdepthdict[bc] = 0
      bcdepthdict[bc] += 1
    inhandle.close()
  
  outfile = open(workpath+'summary_RNA/readcount_barcode_'+sample+'.txt','w') 
  #barcode read_counts UMI_counts
  UMIfile = open(workpath+'summary_RNA/readcount_UMI_'+sample+'.txt','w') 
  #UMI read_count best_barcode best_barcode_freq
  bcdict = {}
  bccount = 0
  for UMI in countdict:
    bccount += 1
    if bccount % 1000000 == 0: print('writing UMI '+str(bccount)+'\ntime: '+str(time.time()-start_time))
    record = Counter(countdict[UMI]).most_common(1)
    bc = record[0][0]
    count = record[0][1]
    depth = len(countdict[UMI])
    UMIfile.write(UMI+'\t'+str(depth)+'\t'+bc+'\t'+str(count/float(depth))+'\n')
    if bc not in bcdict:
      bcdict[bc] = 0
    bcdict[bc] += 1
  for bc in bcdepthdict:
    outfile.write(bc+'\t'+str(bcdepthdict[bc]))
    if bc not in bcdict:
      outfile.write('\t0\n')
    else:
      outfile.write('\t'+str(bcdict[bc])+'\n')
  outfile.close()
  UMIfile.close()
        
if __name__ == '__main__':
  main()
