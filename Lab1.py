from typing import Sequence
from Bio import SeqIO
from Bio.Seq import Seq,reverse_complement, transcribe, back_transcribe, translate
from Bio.SeqRecord import SeqRecord
import numpy as np
import itertools


def findORF(seq,x):
     ORF = []
     for letter in range(x, len(seq), 3):
          endCodon = seq[letter:letter+3]
          if(endCodon == "TAA" or endCodon == "TAG" or endCodon == "TGA" ):
                    for letter1 in range(letter + 3, len(seq),3):
                         StopCodon = seq[letter1:letter1+3]
                         if(StopCodon == "TAA" or StopCodon == "TAG" or StopCodon == "TGA"): 
                              coding = seq[letter+3:letter1+3]
                              letter = letter1
                              if(len(coding)> 100):
                                   ORF.append(coding)
                              break
     return ORF

def findSeq(ORFS,x):
     CodingList = []
     for ORF in ORFS:
          for StartCodon in range(x, len(ORF), 3):
               cod = ORF[StartCodon:StartCodon+3]
               if(cod == "ATG"):
                    coding = ORF[StartCodon:len(ORF)]
                    #print(coding)
                    if(len(coding) > 100):
                         CodingList.append(coding)   
                    break
     return CodingList

def CountCodonFreq(seq):
     totalCodons = len(seq)/3
     Letters = ['A','T','C','G']
     letterProducts = itertools.product(Letters, repeat=3)
     CodonList = []
     for prod in letterProducts:
          CodonList.append(''.join(prod))

     CodonDict = dict.fromkeys(CodonList, 0.0)
     for letter in range(0, len(seq), 3):
          Codon = seq[letter: letter + 3]
          CodonDict[Codon] += 1
     for Codon in CodonDict:
          CodonDict[Codon] = CodonDict[Codon] / totalCodons * 100
     return CodonDict

def CountDicodonFreq(seq):
     if(len(seq) % 6 != 0):
          seq = seq[:-3]
     totalDicodons = len(seq)/6
     Letters = ['A','T','C','G']
     letterProducts = itertools.product(Letters, repeat=6)
     DicodonList = []
     for prod in letterProducts:
          DicodonList.append(''.join(prod))
     
     DicodonDict = dict.fromkeys(DicodonList, 0.0)
     for letter in range(0, len(seq), 6):
          Dicodon = seq[letter: letter + 6]
          DicodonDict[Dicodon] += 1
     for Dicodon in DicodonDict:
          DicodonDict[Dicodon] = DicodonDict[Dicodon] / totalDicodons * 100
     return DicodonDict

def start_stop(seq,x):
     CodingList = []
     for letter in range(x, len(seq),3):
          startCodon = seq[letter:letter+3]
          if(startCodon == "ATG"):
               for letter1 in range(letter, len(seq),3):
                    endCodon = seq[letter1:letter1+3]
                    if(endCodon == "TAA" or endCodon == "TAG" or endCodon == "TGA"):
                         coding = seq[letter:letter1+3]
                         letter = letter1 + 3 
                         if(len(coding) > 100):
                              CodingList.append(coding)
                         break
     return CodingList

def findCoding(Fasta):
     CodingList1 = findSeq(findORF(Fasta.seq,0), 0)
     CodingList2 = findSeq(findORF(Fasta.seq,1), 0)
     CodingList3 = findSeq(findORF(Fasta.seq,2), 0)
     reverse_fasta = Seq(Fasta.seq.reverse_complement())
     CodingList4 = findSeq(findORF(reverse_fasta,0), 0)
     CodingList5 = findSeq(findORF(reverse_fasta,1), 0)
     CodingList6 = findSeq(findORF(reverse_fasta,2), 0)

     Coding1 = Seq('').join(CodingList1)
     Coding2 = Seq('').join(CodingList2)
     Coding3 = Seq('').join(CodingList3)
     Coding4 = Seq('').join(CodingList4)
     Coding5 = Seq('').join(CodingList5)
     Coding6 = Seq('').join(CodingList6)
     AllCoding = Coding1 + Coding2 + Coding3 + Coding4 + Coding5 + Coding6
     return AllCoding

Lactococcus_phage = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\bacterial1.fasta","fasta")
Escherichia_phage = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\bacterial2.fasta","fasta")
Streptococcus_phage = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\bacterial3.fasta","fasta")
Cellulophaga_phage = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\bacterial4.fasta","fasta")
coronavirus = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\mamalian1.fasta","fasta")
adenovirus = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\mamalian2.fasta","fasta")
Variola_virus = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\mamalian3.fasta","fasta")
herpesvirus = SeqIO.read("C:\\Users\\Paulius\\Documents\\Bioinformatics\\viruses\\data\\mamalian4.fasta","fasta")

Codon = CountCodonFreq(findCoding(Lactococcus_phage))

CodonFreq = []
CodonFreq.append(CountCodonFreq(findCoding(Lactococcus_phage)))
CodonFreq.append(CountCodonFreq(findCoding(Escherichia_phage)))
CodonFreq.append(CountCodonFreq(findCoding(Streptococcus_phage)))
CodonFreq.append(CountCodonFreq(findCoding(Cellulophaga_phage)))
CodonFreq.append(CountCodonFreq(findCoding(coronavirus)))
CodonFreq.append(CountCodonFreq(findCoding(adenovirus)))
CodonFreq.append(CountCodonFreq(findCoding(Variola_virus)))
CodonFreq.append(CountCodonFreq(findCoding(herpesvirus)))


Freq = []

CodonMatrix = [[0 for i in range(8)] for j in range (8)]
for i in range(0,8):
     for j in range(0,8):
          for key in Codon:
               Freq.append((CodonFreq[i][key] - CodonFreq[j][key])**2)
          CodonMatrix[i][j] = sum(Freq)
          Freq = []

np.savetxt("CodonFreq.txt", CodonMatrix)

Dicodon = CountDicodonFreq(findCoding(Lactococcus_phage))
DicodonFreq = []
DicodonFreq.append(CountDicodonFreq(findCoding(Lactococcus_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(Escherichia_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(Streptococcus_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(Cellulophaga_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(coronavirus)))
DicodonFreq.append(CountDicodonFreq(findCoding(adenovirus)))
DicodonFreq.append(CountDicodonFreq(findCoding(Variola_virus)))
DicodonFreq.append(CountDicodonFreq(findCoding(herpesvirus)))

Freq = []

DicodonMatrix = [[0 for i in range(8)] for j in range (8)]
for i in range(0,8):
     for j in range(0,8):
          for key in Dicodon:
               Freq.append((DicodonFreq[i][key] - DicodonFreq[j][key])**2)
          DicodonMatrix[i][j] = sum(Freq)
          Freq = []

np.savetxt("DicodonFreq.txt", DicodonMatrix)
