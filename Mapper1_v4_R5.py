#!/usr/bin/env python
#ISS using v4 primers
#reverse read: LTR-genome-Ada-UMI-Ada-Env...
#forward read: Nef-barcode-Env-Ada-UMI-Ada-genome
#goal: truncate viral genome, annotate barcode and UMI in fastq sequence name, write trimmed fastq.
#filters: incompelete reads, missing constant regions, swapped constant regions.
from Bio import SeqIO
from distance import hamming
import sys
from time import time

def main():
  workpath = '/u/scratch/t/tianhao/NovaSeq100622/'
  #070121 did not work well. It is greatly contaminated by BCAmp NL43.
  infile1 = workpath+'split/'+sys.argv[1]
  #test
  #infile1 = workpath+'split/B40-IS12-501712_S12_L002_R1_003.fastq'
  #infile1 = workpath+'split/IS-C03-S25-505701_S31_L001_R1_001.fastq' #020222
  #infile1 = workpath+'split/
  infile2 = infile1.replace('_R1_','_R2_')
  #samplen = infile1.rsplit('_')[-1].rsplit('.')[0]
  samplen = infile1.rsplit('/')[-1].rsplit('_')[-1].rsplit('.')[0]
  #group   = 'R5_'+infile1.rsplit('C03-')[-1].rsplit('-')[0] #020222
  #group   = 'undetermined'
  group   = 'NFNSX_'+infile1.rsplit('/')[-1].rsplit('_')[0]
  outfile1 = open(workpath+'ISS_genome/'+group+'_R1_'+samplen+'.fastq','w')
  outfile2 = open(workpath+'ISS_genome/'+group+'_R2_'+samplen+'.fastq','w')
  unintefile = open(workpath+'ISS_unintegrated/'+group+'_'+samplen+'.fastq','w')
  plasmidfile = open(workpath+'ISS_plasmid/'+group+'_plasmid_'+samplen+'.txt','w')
  inhandle1 = open(infile1); handle1 = SeqIO.parse(inhandle1,'fastq')
  inhandle2 = open(infile2); handle2 = SeqIO.parse(inhandle2,'fastq')
  #read2(R5): TAGTCAGTGTGGAAAATCTCTAGCA: genome: ATTGAGGTTTGCAGTTG(Ada) 
  #read1(R5): TTTTGACCACTTGCCACCCAT(Nef5):21ANN barcode:CAAAGCCCTTTCCAAGCC(Env3):GAATTC(EcoRI):10N UMI:CAACTGCAAACCTCAAT(FL):genome:TGCTAGAGATTTTCCAC(LTRRC)
  logfile = open(workpath+'ISS_log/adapter_mapping_'+group+'_'+samplen+'.txt','w')
  count = [0]*8; start_time = time()
  for record1 in handle1:
    if count[0] % 1000000 == 0:
      print('Time: '+str(time()-start_time)+'s\tCount table:'+str(count))
      #10s for 10000 records. For 0.63M record, it takes ~10min. request 3 hours for each 10 samples
    count[0] += 1
    record2 = next(handle2)
    pltr = 'GTCAGTGTGGAAAATCTCT'; padaf = 'ATTGAGGTTTGCAGTTG'
    pnef = 'TTGACCACTTGCCACCCAT'; penv = 'CAAAGCCCTTTCCAAGCCCT'; peco = 'GAATTC'; pada = 'CAACTGCAAACCTCAAT'; pltrc = 'AGAGATTTTCCACAC'
    rpla = 'CCCAGGAGGTAGAGGTTGCAGTGAGCCAA'
    seq1 = str(record1.seq); seq2 = str(record2.seq)
    f_offset,f_offset2 = map_read1(seq2,pltr,padaf)
    r_offset,r_offset2,r_offset3,r_offset4,r_offset5,r_offset6 = map_read2(seq1,pnef,penv,peco,pada,pltrc)
    if f_offset == 'NA' or r_offset == 'NA' or r_offset3 == 'NA': 
      count[2] += 1
      #print(f_offset,f_offset2)
      #print(r_offset,r_offset2,r_offset3,r_offset4,r_offset5,r_offset6)
      continue
    barcode = seq1[r_offset:r_offset2]
    UMI     = seq1[r_offset3:r_offset4]
    if len(barcode) > 60 or len(UMI) > 40 or len(barcode) == 0 or len(UMI) == 0: 
      count[5] += 1
      logfile.write(record1.format('fastq')+record2.format('fastq')+'_'.join([str(x) for x in [f_offset,f_offset2,r_offset,r_offset2,r_offset3,r_offset4,r_offset5,r_offset6]])+'\n')
      continue
    if f_offset2 != 'NA':
      record2 = record2[f_offset:f_offset2]
    else:
      record2 = record2[f_offset:]
    #print(record1.id)
    #if '1138:24840:9878' in str(record1.id):
    #  print(f_offset,f_offset2)
    #  print(r_offset,r_offset2,r_offset3,r_offset4,r_offset5,r_offset6)
    #  break
    if r_offset6 != 'NA':
      record1 = record1[r_offset5:r_offset6]
    else:
      record1 = record1[r_offset5:]
    if len(record2.seq) < 10 or len(record1.seq) < 10:
      unintefile.write(UMI+'\t'+barcode+'\t'+str(record2.seq)+'\n')
      count[6] += 1
      continue
    seq2 = str(record2.seq)
    _tmp1,_tmp2 = map_oligo(seq2,rpla)
    if _tmp2 > 0:
      plasmidfile.write(UMI+'\t'+barcode+'\t'+str(record2.seq)+'\n')
      count[7] += 1
      continue
    record1.id += (':'+UMI+':'+barcode)
    record2.id += (':'+UMI+':'+barcode)
    record1.description = ''; record2.description = ''
    count[1] += 1
    outfile1.write(record1.format('fastq'))
    outfile2.write(record2.format('fastq'))
  logfile.write('total read:'+str(count[0])+'\n')
  logfile.write('effective read:'+str(count[1])+'\n')
  logfile.write('problematic constant region:'+str(count[2])+'\n')
  logfile.write('SapI read:'+str(count[3])+'\n') 
  logfile.write('EcoRI read:'+str(count[4])+'\n') 
  logfile.write('wierd mapping read:'+str(count[5])+'\n') 
  logfile.write('short read:'+str(count[6])+'\n')
  logfile.write('plasmid read:'+str(count[7])+'\n')
  inhandle1.close(); inhandle2.close(); outfile1.close(); outfile2.close();logfile.close();unintefile.close();plasmidfile.close()
  
def map_read1(seq,pltr,padaf):
  _tmp, f_offset = map_oligo(seq,pltr)
  f_offset2,_tmp = map_oligo(seq[f_offset:],padaf)
  if f_offset == -1: 
    return 'NA','NA'
  if _tmp == -1: 
    f_offset2 = 'NA'
  else:
    f_offset2 += f_offset
  return f_offset,f_offset2

def map_read2(seq,pnef,penv,peco,pada,pltrc):
  _tmp, r_offset = map_oligo(seq,pnef)
  if r_offset == -1: 
    #print('pnef')
    return 'NA','NA','NA','NA','NA','NA'
  r_offset2, r_offset3 = map_oligo(seq[r_offset:],penv+peco)
  if r_offset2 == -1 or r_offset3 == -1:
    #print('penv+peco')
    return 'NA','NA','NA','NA','NA','NA'
  r_offset2 += r_offset
  r_offset3 += r_offset
  r_offset4, r_offset5 = map_oligo(seq[r_offset3:],pada)
  if r_offset4 == -1 or r_offset5 == -1:
    #print('pada')
    return 'NA','NA','NA','NA','NA','NA'
  r_offset4 += r_offset3
  r_offset5 += r_offset3
  r_offset6, _tmp = map_oligo(seq[r_offset5:],pltrc)
  if r_offset6 == -1:
    r_offset6 = 'NA'
  else:
    r_offset6 += r_offset5
  return r_offset,r_offset2,r_offset3,r_offset4,r_offset5,r_offset6

def map_oligo(seq,refseq):
  _offset = [-1,-1]
  for offset in range(250):
    for ref_trim in range(2):
      tref = refseq[ref_trim:]
      tarseq = seq[offset:offset+len(tref)]
      if len(tarseq) != len(tref): break
      if hamming(tarseq,tref) < 4:
        _offset[1] = offset+len(tref)
        _offset[0] = offset
        break
    if _offset[0] != -1:
      break
  return _offset[0],_offset[1]

if __name__ == '__main__':
  main()
