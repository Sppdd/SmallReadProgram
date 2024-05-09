#!/usr/bin/python 
#python program location

dna = input('Enter your DNA Squance, please:')

no_c = dna.count('c')
no_g = dna.count('g')

dna_length = len(dna)

gc_percent = (no_c + no_g) * 100.0 / dna_length

print(gc_percent, '%')

def hasStopCodon(dna,frame = 0):
  CodonFound = False
  StopCodons = ['tag', 'tag','taa']
  for i in range(frame,len(dna),3):
   codon=dna [i:i+3].lower()
   if codon in StopCodons :
      CodonFound = True
      break

  return CodonFound