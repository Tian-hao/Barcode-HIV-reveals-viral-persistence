#!/usr/bin/env python
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from time import time
from distance import hamming
#1hour per 1M records (4M lines)
#read1: BCF(GGCTTGGAAAGGATTTTGCTATAA)-NNCNNCNNCNNCNNCNNCNNC-rcBCR(TATAAGATGGGTGGCAAGTG)-21N UMI-rchd(CTCGAGGATCTAGGCACTTG)
#read2: handle(CAAGTGCCTAGATCCTCGAG)-21N UMI-BCR(CACTTGCCACCCATCTTATA)-GNNGNNGNNGNNGNNGNNGNN-rcBCF(TTATAGCAAAATCCTTTCCAAGCC)
workpath = '/u/scratch/t/tianhao/NovaSeq040122/'
priseqs = {'BCF':'GGCTTGGAAAGGATTTTGCTATAA',
'rcBCR':'TATAAGATGGGTGGCAAGTG',
'rchd':'CTCGAGGATCTAGGCACTTG',
'handle':'CAAGTGCCTAGATCCTCGAG',
'BCR':'CACTTGCCACCCATCTTATA',
'rcBCF':'TTATAGCAAAATCCTTTCCAAGCC'} 
Flist  = ['BCF','rcBCR','rchd']
Rlist  = ['handle','BCR','rcBCF']

def main(): 
  infile1 = workpath+'split/'+sys.argv[1]
  infile2 = infile1.replace('_R1_','_R2_')
  sample = infile1.rsplit('/')[-1].rsplit('.fastq')[0].rsplit('_')
  sample = sample[0]+'_'+sample[-1]
  outfile = open(workpath+'ISS_RNA/RNA_barcode_'+sample+'.txt','w')
  inhandle1 = open(infile1); inhandle2 = open(infile2)
  handle1 = SeqIO.parse(inhandle1,'fastq'); handle2 = SeqIO.parse(inhandle2,'fastq')
  errcount = [0,0,0,0,0,0]; readcount = 0
  #0: unpaired reads
  #1. unmapped
  #2. only 1 primer mapped
  #3. primer order wrong
  #4. low quality
  #5. good count
  start_time = time()
  for record1 in handle1:
    record2 = next(handle2)
    readcount += 1
    if readcount % 1000000 == 0:
      print(errcount,readcount,'Time: '+str(time()-start_time))
    lib = sample
    #map barcode
    seq1 = str(record1.seq)
    seq2 = str(record2.seq)
    direction1 = mapdir(seq1)
    direction2 = mapdir(seq2)
    if direction1 == 'NA' or direction2 == 'NA': errcount[1] += 1; continue
    if direction1 == direction2: errcount[3] += 1; continue
    offsets1 = mapbc(seq1,direction1)
    offsets2 = mapbc(seq2,direction2)
    if len(offsets1) == 0 or len(offsets2) == 0: errcount[2] += 1; continue
    if offsets1[2] < offsets1[1] or offsets1[1] < offsets1[0] or offsets2[2] < offsets2[1] or offsets2[1] < offsets2[0]: errcount[3] += 1; continue
    if direction1 == 'F':
      bcseq1,bcqual1,UMIseq1,UMIqual1 = extractseq(record1,offsets1,Flist)
      bcseq2,bcqual2,UMIseq2,UMIqual2 = extractseq(record2,offsets2,Rlist)
    if direction1 == 'R':
      bcseq2,bcqual2,UMIseq2,UMIqual2 = extractseq(record1,offsets1,Rlist)
      bcseq1,bcqual1,UMIseq1,UMIqual1 = extractseq(record2,offsets2,Flist)      
    if len(bcseq1) != len(bcseq2) or 'N' in bcseq1 or 'N' in bcseq2: errcount[4] += 1; continue
    #write barcode
    low_qs_flag, unpair_flag, bcseq = consense2(bcseq1,bcseq2,bcqual1,bcqual2)
    if low_qs_flag == 1: errcount[4] += 1; continue
    if unpair_flag == 1: errcount[0] += 1; continue
    if len(UMIseq1) > len(UMIseq2): 
      UMIseq1 = UMIseq1[:len(UMIseq2)]
      UMIqual1 = UMIqual1[:len(UMIseq2)]
    if len(UMIseq1) < len(UMIseq2): 
      UMI = str(UMIseq2)
    else:
      low_qs_flag, unpair_flag, UMI = consense2(UMIseq1,UMIseq2,UMIqual1,UMIqual2)
    if low_qs_flag == 1: errcount[4] += 1; UMI = 'NA'
    if unpair_flag == 1: errcount[0] += 1; UMI = 'NA'
    outfile.write(UMI+'\t'+bcseq+'\n')
    errcount[5] += 1
  inhandle1.close()
  inhandle2.close()
  outfile.write('Unpaired reads: '+str(errcount[0])+'\n')
  outfile.write('Unmapped reads: '+str(errcount[1])+'\n')
  outfile.write('1-Primer reads: '+str(errcount[2])+'\n')
  outfile.write('Primers not in order reads: '+str(errcount[3])+'\n')
  outfile.write('Low quanlity reads: '+str(errcount[4])+'\n')
  outfile.close()

def consense2(bcseq1,bcseq2,bcqual1,bcqual2):
  low_qs_flag = 0; unpair_flag = 0
  bcseq = ''
  for n in range(len(bcseq1)):
    if bcseq1[n] == bcseq2[n]: bcseq += bcseq1[n]
    elif bcseq1[n] != bcseq2[n] and bcqual1[n] >= 30 and bcqual2[n] < 30: bcseq += bcseq1[n]
    elif bcseq1[n] != bcseq2[n] and bcqual1[n] < 30 and bcqual2[n] >= 30: bcseq += bcseq2[n]
    elif bcseq1[n] != bcseq2[n] and bcqual1[n] >= 30 and bcqual2[n] >= 30: unpair_flag = 1; break
    elif bcseq1[n] != bcseq2[n] and bcqual1[n] < 30 and bcqual2[n] < 30: low_qs_flag = 1; break
  return low_qs_flag, unpair_flag, bcseq

def extractseq(record,offsets,prilist):
  bcseq   = record[offsets[0]+len(priseqs[prilist[0]]):offsets[1]].seq
  bcqual  = record[offsets[0]+len(priseqs[prilist[0]]):offsets[1]].letter_annotations["phred_quality"]
  UMIseq  = record[offsets[1]+len(priseqs[prilist[1]]):offsets[2]].seq
  UMIqual = record[offsets[1]+len(priseqs[prilist[1]]):offsets[2]].letter_annotations["phred_quality"]
  if 'handle' in prilist:
    t       = bcqual[::-1]
    bcqual  = UMIqual[::-1]
    UMIqual = t
    t       = bcseq.reverse_complement()
    bcseq   = UMIseq.reverse_complement()
    UMIseq  = t
  return bcseq,bcqual,UMIseq,UMIqual
    
def mapdir(seq):
  for direction in ['BCF','handle']:
    priseq = priseqs[direction]
    offset = mapseq(seq,priseq)
    if offset != 'NA':
      if direction == 'BCF': return 'F'
      if direction == 'handle': return 'R' 
  return 'NA'

def mapbc(seq,direction):
  offsets = []
  if direction == 'F':
    #seq = seq[20:]
    for pri in Flist:
      priseq = priseqs[pri]
      offset = mapseq(seq,priseq)
      if offset == 'NA': return []
      offsets.append(offset)
  if direction == 'R':
    for pri in Rlist:
      priseq = priseqs[pri]
      offset = mapseq(seq,priseq)
      if offset == 'NA': return []
      offsets.append(offset)
  return offsets

def mapseq(seq,pri):
  for offset in range(150):
    qseq = seq[offset:offset+len(pri)]
    if len(qseq) < len(pri): break
    if hamming(qseq,pri) <= 3:
      return offset
  return 'NA'

if __name__ == '__main__':
  main()
