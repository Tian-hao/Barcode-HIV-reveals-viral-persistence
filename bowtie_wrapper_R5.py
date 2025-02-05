#!/usr/bin/env python
import os
import glob
import sys

def main():
  workpath = '/u/scratch/t/tianhao/NovaSeq100622/'
  homepath = '/u/home/t/tianhao/ISS/'
  infiles = sorted(glob.glob(workpath+'ISS_genome/NFNSX_*_R1_*.fastq'))
  samplelist = [x.rsplit('/')[-1].rsplit('_')[0]+'_'+x.rsplit('/')[-1].rsplit('_')[1] for x in infiles]
  samplelist = list(set(samplelist))
  print(len(samplelist))
  for sample in samplelist:
    os.system('cat '+workpath+'ISS_genome/'+sample+'_R1_*.fastq > '+workpath+'ISS_genome/combined_'+sample+'_R1.fastq')
    os.system('cat '+workpath+'ISS_genome/'+sample+'_R2_*.fastq > '+workpath+'ISS_genome/combined_'+sample+'_R2.fastq')
    infile1 = workpath+'ISS_genome/combined_'+sample+'_R1.fastq'
    infile2 = infile1.replace('_R1','_R2')
    for j in ['hg38','HIV']: 
      outfile = workpath+'ISS_bowtie/'+j+'/'+sample+'.sam'
      metfile = workpath+'ISS_log/bowtie_mapping_'+j+'_'+sample+'.txt'
      bashfile = open(homepath+'/script/shell_tmp/bowtie_wrapper_'+j+'_'+sample+'.sh','w')
      bashfile.write('#!/bin/bash\n')
      if j == 'hg38': reffile = '/u/scratch/t/tianhao/NovaSeq100622/ref/hg38'
      #if j == 'HIV': reffile = '~/ISS/ref/HIV'
      if j == 'HIV': reffile = '/u/scratch/t/tianhao/NovaSeq100622/ref/NFNSX'
      bashfile.write('bowtie2 --very-sensitive -x '+reffile+' -1 '+infile1+' -2 '+infile2+' -S '+outfile+' --met-file '+metfile)
      #bashfile.write('bowtie2 --very-sensitive -x '+reffile+' -U '+infile1+' -S '+outfile)
      bashfile.close()
      os.system('chmod 777 '+homepath+'script/shell_tmp/bowtie_wrapper_'+j+'_'+sample+'.sh')
      os.system('qsub -cwd -V -N PJbowtie -l h_data=6144M,h_rt=3:00:00 '+homepath+'/script/shell_tmp/bowtie_wrapper_'+j+'_'+sample+'.sh')

if __name__ == '__main__':
  main()
