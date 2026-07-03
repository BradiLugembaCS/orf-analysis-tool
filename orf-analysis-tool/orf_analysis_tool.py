#C22348001 - Bioinformatics Assignment

"""
3. Using a suitable example such as the e. coli pal gene sequence. Write a
python program(s) that 1:
a. Show the frame number; the nucleotide and amino acid start and
stop positions; the length of the ORFs in nucleotide/ amino
positions and display the DNA/AA sequence for all potential
Realistic O.R.F.

"""

# used to open windows explorer to read files
import tkinter as tk
from tkinter import filedialog
import sys


# DNA <-> AA translation table: CodonTable
CodonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}


#************* Read the contents of the fasta file *********************

def FileRead():

    """
        this function opens a fasta file
        it reads the descriptor line and all DNA lines
        removes '\n' and returns descriptor line and DNA sequence
    """

    Data = ""
    DesLine = ""
    FileContents = []

    print("creating file explorer window minimise spyder to see it....")
    root = tk.Tk()
    root.wm_withdraw()
    FileName = filedialog.askopenfilename(filetypes=[('All files', '*.*')])
    root.destroy()

    try:
        Fp1 = open(FileName, 'r')
        DesLine = Fp1.readline().rstrip('\n')
        Data = Fp1.read()
        Fp1.close()

    except IOError:
        print("error unable to read file or file does not exist!!!")
        print("Exiting the program")
        input("press return")
        sys.exit(1)

    # remove end of line characters and spaces
    ListSeq = Data.split('\n')
    DnaSeq = ('').join(ListSeq)
    DnaSeq = DnaSeq.replace(" ", "").upper()

    FileContents.append(DesLine)
    FileContents.append(DnaSeq)

    return FileContents


#*********************************** Translate a sequence from a given start position ***************************

def Translate(DnaSequence, RFNumber):

    """
        this function translates a DNA sequence into an amino acid sequence
        starting at the reading frame number given
    """

    AminoAcidSeq = ''

    for n in range(RFNumber, len(DnaSequence), 3):
        codon = DnaSequence[n:n+3]

        if len(codon) == 3 and codon in CodonTable:
            AminoAcid = CodonTable[codon]
            AminoAcidSeq += AminoAcid

    return AminoAcidSeq


# ********************** function to get the compliment of a DNA sequence *********************************

def Compliment(DnaSeq):

    ComplimentSeq = ''

    for index in range(0, len(DnaSeq)):
        if DnaSeq[index] == 'T':
            ComplimentSeq += 'A'
        if DnaSeq[index] == 'A':
            ComplimentSeq += 'T'
        if DnaSeq[index] == 'C':
            ComplimentSeq += 'G'
        if DnaSeq[index] == 'G':
            ComplimentSeq += 'C'

    return ComplimentSeq


#******************** Find ORFs in one amino acid strand *******************

def AllORF(AAStrand, DNASequenceRF, FrameLabel, RFOffset, MinAA):

    """
        this method finds all potential ORF M followed by *
        it stores important ORF information in a list

        it takes:
            amino acid strand
            DNA strand for that same reading frame
            frame label
            frame offset
            minimum AA cutoff

        it returns a list of ORFs
    """

    ORFList = []

    index = 0
    while index < len(AAStrand):

        if AAStrand[index] == 'M':         # found start codon in AA sequence
            startAA = index + 1            # amino acid positions shown as 1-based
            ORF = ""
            innerIndex = index

            # collect amino acids until stop codon or end of strand
            while innerIndex < len(AAStrand) and AAStrand[innerIndex] != '*':
                ORF += AAStrand[innerIndex]
                innerIndex += 1

            # only keep if a stop codon was found
            if innerIndex < len(AAStrand) and AAStrand[innerIndex] == '*':
                stopAA = innerIndex + 1    # stop codon position in AA strand
                ORFLengthAA = len(ORF)
                ORFLengthNT = ORFLengthAA * 3

                # realistic ORF cutoff
                if ORFLengthAA >= MinAA:

                    # nucleotide positions on the reading frame DNA sequence
                    startNT = (startAA - 1) * 3 + 1
                    stopNT = startNT + ORFLengthNT - 1

                    DNAORF = DNASequenceRF[startNT - 1: stopNT]

                    ORFInfo = {
                        "frame": FrameLabel,
                        "aa_start": startAA,
                        "aa_stop": stopAA,
                        "aa_length": ORFLengthAA,
                        "nt_start": startNT + RFOffset,
                        "nt_stop": stopNT + RFOffset,
                        "nt_length": ORFLengthNT,
                        "dna_seq": DNAORF,
                        "aa_seq": ORF
                    }

                    ORFList.append(ORFInfo)

            # continue looking from the amino acid after the current start
            index += 1

        else:
            index += 1

    return ORFList


#**************************** print ORFs nicely **************************************

def PrintORFs(ORFList):

    if len(ORFList) == 0:
        print("\nNo realistic ORFs found in this reading frame.")
        return

    for count, ORF in enumerate(ORFList, start=1):

        print("\n==============================================================")
        print("ORF Number: {:d}".format(count))
        print("Frame Number: {:s}".format(ORF["frame"]))
        print("Nucleotide Start Position: {:d}".format(ORF["nt_start"]))
        print("Nucleotide Stop Position : {:d}".format(ORF["nt_stop"]))
        print("Amino Acid Start Position: {:d}".format(ORF["aa_start"]))
        print("Amino Acid Stop Position : {:d}".format(ORF["aa_stop"]))
        print("ORF Length in Nucleotides: {:d}".format(ORF["nt_length"]))
        print("ORF Length in Amino Acids: {:d}".format(ORF["aa_length"]))
        print("\nDNA Sequence:")
        print(ORF["dna_seq"])
        print("\nAmino Acid Sequence:")
        print(ORF["aa_seq"])
        print("==============================================================")


#*******************************  the driver or main function ********************************************************

def main():

    DesLine = ''
    DnaSeq = ''
    FileContents = []

    # read fasta file
    FileContents = FileRead()
    DesLine = FileContents[0]
    DnaSeq = FileContents[1]


    print("\nDescriptor line is:")
    print(DesLine)

    print("\nThe DNA sequence is:")
    print(DnaSeq)

    # ask user for minimum realistic ORF size
    try:
        MinAA = int(input("\nEnter the minimum ORF size in amino acids (e.g. 30): "))
    except ValueError:
        print("invalid entry... using default value of 20")
        MinAA = 20

    # reverse compliment
    ComplimentSeq = Compliment(DnaSeq)
    ReverseComplimentSeq = ComplimentSeq[::-1]

    print("\n*************************  All potential realistic ORFs of the PRIMARY STRAND  *************************\n")

    for RFNumber in range(0, 3):
        FrameLabel = "+{:d}".format(RFNumber + 1)
        DnaSequenceRF = DnaSeq[RFNumber:]
        AASequence = Translate(DnaSeq, RFNumber)

        print("\n************************************* reading frame {:s} *******************************".format(FrameLabel))
        print("The amino acid sequence is:\n")
        print(AASequence)

        ORFList = AllORF(AASequence, DnaSequenceRF, FrameLabel, RFNumber, MinAA)
        PrintORFs(ORFList)

    print("\n*************************  All potential realistic ORFs of the REVERSE COMPLIMENTARY STRAND  *************************\n")

    for RFNumber in range(0, 3):
        FrameLabel = "-{:d}".format(RFNumber + 1)
        DnaSequenceRF = ReverseComplimentSeq[RFNumber:]
        AASequence = Translate(ReverseComplimentSeq, RFNumber)

        print("\n************************************* reading frame {:s} *******************************".format(FrameLabel))
        print("The amino acid sequence is:\n")
        print(AASequence)

        ORFList = AllORF(AASequence, DnaSequenceRF, FrameLabel, RFNumber, MinAA)
        PrintORFs(ORFList)

    input("\nPress return to finish....")


#**************** execute program **************************

main()