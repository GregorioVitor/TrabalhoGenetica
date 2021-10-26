import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("-F", "--fastaFile", required=False, help= "Path of the fasta file")
ap.add_argument("-S", "--sequence", required=False, help= "Sequence in plain text")
args = vars(ap.parse_args())


def achar_orf_F(seq, pos, j):
    ORF = seq.find("ATG",pos)
    frame=ORF+1
    while ORF != -1:
        for i in range(ORF, len(seq), 3):
            if(seq[i:i+3] == 'TAA' or seq[i:i+3] == 'TAG' or seq[i:i+3] == 'TGA'):
                if(ORF%3==0):
                    frame=1
                elif(ORF%3==1):
                    frame=2
                else:
                    frame=3
                print("ORF: "+ str(j)+ " - Frame: "+ str(frame))
                print("ORF está na posição "+str(ORF+1)+ "; "+str(ORF+2)+ " e "+str(ORF+3))
                print("Codon de parada está na posição "+str(i+1)+"; "+str(i+2)+" e "+str(i+3))
                print("Sequencia: "+ str(seq[ORF:i+3]))
                j+=1
                print()
                achar_orf_F(seq, i, j)
                ORF=-1
                break
            elif (i >= len(seq)-5):
                achar_orf_F(seq, ORF+1, j)
                ORF=-1
                break
            elif (i >= len(seq)-3):
                ORF=-1
                break  
            
if(args['fastaFile']!=None):
    i=0
    for seq_record in SeqIO.parse(args['fastaFile'], "fasta"):
        print("Sequencia "+ str(i))
        achar_orf_F(seq_record.seq, 0, 1)
        print()
        i+=1
elif(args['sequence']!=None):
    achar_orf_F(args['sequence'], 0, 1)
else:
    print("Invalid Syntax, use -F for fasta Files, or -S for sequence in plain text")    

