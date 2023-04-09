import argparse
#pKa data from EMBOSS
pKn={"C":8.5,"D":3.9,"E":4.1,"Y":10.1}
pKp={"R":12.5,"K":10.8,"H":6.5}


parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create\
    the pH-dependent charge file from fasta file. Written by Xinyu Gu."
)

parser.add_argument("-f", "--fasta", help="fasta file")
parser.add_argument("--pH", help="pH value")
parser.add_argument("-c", "--charge", help="charge file")
args = parser.parse_args()


fastaFile=args.fasta
chargeFile=args.charge
pH=float(args.pH)


seq = ""
with open(fastaFile) as f:
    for line in f:
        if line[0] == ">":
            pass
        elif(line == "\n"):
            pass
        else:
            # print(line)
            seq += line.strip("\n")



with open(chargeFile, "w") as out:
    for i in seq:
        if i=="C" or i=="D" or i=="E" or i=="Y":
           charge=-1/(1+10**(pKn[i]-pH))
           out.write(f"{charge:.2f}\n")
        elif i=="R" or i=="K" or i=="H":
           charge=1/(1+10**(pH-pKp[i]))
           out.write(f"{charge:.2f}\n")
        else:
           out.write("0\n")

out.close()
