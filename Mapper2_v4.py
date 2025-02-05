#!/usr/bin/env python
#read sams (human and HIV), plasmids and unintegrated
#integrate 4 files as archiving data
#write a summation linkage file
#use with mywrappers and use sample name as inputs
import re
import sys
import glob

workpath = '/u/scratch/t/tianhao/NovaSeq100622/'

def main():
  infiles = sorted(glob.glob(workpath+'ISS_plasmid/NFNSX_*.txt'))
  samplelist = []
  for infile in infiles:
    tag = infile.rsplit('/')[-1].rsplit('_')
    tag = tag[0]+'_'+tag[1]
    sample = tag
    if sample not in samplelist:
      print(sample)
      link_record(tag,sample)
      samplelist.append(sample)

def link_record(tag,sample): 
  #tag = sys.argv[1]
  #sample = sys.argv[2]
  archfile  = open(workpath+'ISS_archive/record_'+sample+'.txt','w')
  linkfile  = open(workpath+'ISS_linkage/linkage_'+sample+'.txt','w')
  sumfile   = open(workpath+'ISS_archive/summary_'+sample+'.txt','w')
  plasfiles = sorted(glob.glob(workpath+'ISS_plasmid/'+tag+'_*.txt'))
  uninfiles = sorted(glob.glob(workpath+'ISS_unintegrated/'+tag+'_*.fastq'))
  hg38files = sorted(glob.glob(workpath+'ISS_bowtie/hg38/'+tag+'.sam'))
  nl43files = sorted(glob.glob(workpath+'ISS_bowtie/HIV/'+tag+'.sam'))
  #archive file: record by line
  #UMI	barcode	type	annotation
  #linkage file: UMI by line
  #UMI	UMIcount	barcode	top_barcode_freq	2nd_barcode_freq	type	annotation	top_annotation_freq	2nd_annotation_freq
  #summary file: 
  #mapped to each type, total UMI, total unique barcode
  
  linkdict = {}; logdict = {'unmapped':0,'unpaired':0,'positive_strand':0,'negative_strand':0}
  for infile in hg38files:
    linkdict,logdict = readsam2(linkdict,logdict,infile,'hg38',archfile)
  for infile in nl43files:
    linkdict,logdict = readsam2(linkdict,logdict,infile,'HIV', archfile)
  for infile in uninfiles:
    linkdict,logdict = readtxt(linkdict,logdict,infile,'free',archfile)
  for infile in plasfiles:
    linkdict,logdict = readtxt(linkdict,logdict,infile,'plasmid',archfile)

  bclist = []; annotdict = {}; depth = 0; depdict = {}
  for UMI in linkdict:
    topbc,bctopfreq,bcsecfreq  = find_most([x[0] for x in linkdict[UMI]])
    topannot,atopfreq,asecfreq = find_most([x[1] for x in linkdict[UMI]])
    atype = topannot.rsplit(':')[0]
    topa  = topannot.rsplit(':')[1:]
    topa  = ':'.join(topa)
    linkfile.write(UMI+'\t'+str(len(linkdict[UMI]))+'\t'+topbc+'\t'+str(bctopfreq)+'\t'+str(bcsecfreq)+'\t'+atype+'\t'+topa+'\t'+str(atopfreq)+'\t'+str(asecfreq)+'\n')
    bclist.append(topbc)
    Udep  = len(linkdict[UMI])
    depth += Udep
    if atype not in annotdict: annotdict[atype] = 0
    annotdict[atype] += 1
    if Udep not in depdict: depdict[Udep] = 0
    depdict[Udep] += 1  
  
  sumfile.write('total read: '+str(depth)+'\n')
  sumfile.write('total UMI: '+str(len(linkdict))+'\n')
  sumfile.write('unique barcodes: '+str(len(set(bclist)))+'\n')
  if 'hg38' in annotdict:
    sumfile.write('hg38 UMI: '+str(annotdict['hg38'])+'\n')
  if 'HIV' in annotdict:
    sumfile.write('HIV UMI: '+str(annotdict['HIV'])+'\n')
  if 'free' in annotdict:
    sumfile.write('unintegrated UMI: '+str(annotdict['free'])+'\n')
  if 'plasmid' in annotdict:
    sumfile.write('plasmid UMI: '+str(annotdict['plasmid'])+'\n')
  sumfile.write('unpaired reads:'+str(logdict['unpaired'])+'\n')
  sumfile.write('positive strand mapped reads:'+str(logdict['positive_strand'])+'\n')
  sumfile.write('negative strand mapped reads:'+str(logdict['negative_strand'])+'\n')
  sumfile.write('========UMI depth distribution=======\n')
  sumfile.write('UMI_read_depth\tunique_UMI_count\n')
  for Udep in depdict:
    sumfile.write(str(Udep)+'\t'+str(depdict[Udep])+'\n')
  sumfile.close()

def find_most(bclist):
  bcdict = {}
  for bc in bclist:
    if bc not in bcdict: bcdict[bc] = 0
    bcdict[bc] += 1
  sortbcdict = sorted((v,k) for (k,v) in bcdict.items())
  if len(sortbcdict) == 1:
    return sortbcdict[0][1], 1, 'NA'
  else:
    return sortbcdict[-1][1], float(sortbcdict[-1][0])/len(bclist), float(sortbcdict[-2][0])/len(bclist)

def readtxt(linkdict,logdict,infile,txttype,archfile):
  txtfile = open(infile)
  for line in txtfile:
    line    = line.rstrip('\n').rsplit('\t')
    bc      = line[1]
    UMI     = line[0]
    annot   = line[2]
    archfile.write(UMI+'\t'+bc+'\t'+txttype+'\t'+annot+'\n')
    if UMI not in linkdict:
      linkdict[UMI] = []
    linkdict[UMI].append((bc,txttype+':'+annot))
  txtfile.close()
  return linkdict,logdict

def readsam(linkdict,logdict,infile,samtype,archfile): #legacy
  samfile = open(infile)
  for line in samfile:
    if line[0] == '@': continue
    line    = line.rstrip().rsplit('\t')
    UMI     = line[0].rsplit(':')[-2]
    bc      = line[0].rsplit(':')[-1]
    qual    = line[4]
    chrom   = line[2]
    pos     = line[3]
    flag    = bin(int(line[1]))
    if flag[-3] == '1': continue
    if qual == '1': continue
    if len(flag) >= 5 and flag[-5] == '1':
      strand = '-'
    else:
      strand = '+'
    annot = chrom+':'+pos+':'+strand
    archfile.write(UMI+'\t'+bc+'\t'+samtype+'\t'+annot+'\n')
    if UMI not in linkdict:
      linkdict[UMI] = []
    linkdict[UMI].append((bc,samtype+':'+annot))
  samfile.close()
  return linkdict,logdict

def flag2strand(flag):
    if len(flag) > 5 and flag[-5] == '1':
      return '-'
    else:
      return '+'

def readsam2(linkdict,logdict,infile,samtype,archfile):
  #This is for paired end read 
  lc = 0
  samfile = open(infile)
  for line in samfile:
    if line[0] == '@': continue
    lc += 1
    if lc % 2 == 1: 
      line1 = line
      continue 
    line1   = line1.rstrip().rsplit('\t')
    line2   = line.rstrip().rsplit('\t')
    if line1[0] != line2[0]:
      line1 = '\t'.join(line[2])
      lc -= 1
      continue
    UMI     = line1[0].rsplit(':')[-2]
    bc      = line1[0].rsplit(':')[-1]
    #In the new design, R2 is forward read after LTR.
    flag1   = bin(int(line1[1]))
    flag2   = bin(int(line2[1]))
    if len(flag1) > 7 and flag1[-8] == '1':
      read2 = line1; flag2 = bin(int(line1[1]))
      read1 = line2; flag1 = bin(int(line2[1]))
    else:
      read1 = line1
      read2 = line2 
    qual1   = int(line1[4])
    qual2   = int(line2[4])
    if qual2 >= qual1:
      if flag2[-3] == '1' or qual2 == 1: 
        logdict['unmapped'] += 1
        continue
      strand = flag2strand(flag2)
      if strand == '+':
        pos = read2[3]
      else:
        cigar = read2[5]
        mlength = re.findall(r"(\d+)M", cigar)
        mlength = sum([int(x) for x in mlength])
        pos = str(int(read2[3])+mlength)
    else:
      if flag1[-3] == '1' or qual1 == 1:
        logdict['unmapped'] += 1
        continue
      strand = flag2strand(flag1)
      if strand == '+':
        strand = '-'
        cigar = read1[5]
        mlength = re.findall(r"(\d+)M", cigar)
        mlength = sum([int(x) for x in mlength])
        pos = str(int(read1[3])+mlength)
      else:
        strand = '+'
        pos = read1[3]
    if strand == '+': 
      logdict['positive_strand'] += 1
    else:
      logdict['negative_strand'] += 1
    chrom1 = read1[2]
    chrom2 = read2[2]
    if chrom1 != chrom2 and chrom1 != '*' and chrom2 != '*':
      logdict['unpaired'] += 1
      continue
    else:
      chrom = chrom2
    annot = chrom+':'+pos+':'+strand
    archfile.write(UMI+'\t'+bc+'\t'+samtype+'\t'+annot+'\n')
    if UMI not in linkdict:
      linkdict[UMI] = []
    linkdict[UMI].append((bc,samtype+':'+annot))
  samfile.close()
  return linkdict,logdict

if __name__ == '__main__':
  main()
