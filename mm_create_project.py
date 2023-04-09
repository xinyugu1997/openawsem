#!python
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.insert(0, OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein, \
            do: python3 ~/OPENAWSEM_LOCATION/mm_create_project.py")
parser.add_argument("-c", "--chain", default="-1", help="chains to be simulated, could be for example 'abc'.")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--frag", action="store_true", default=False)
parser.add_argument("--extended", action="store_true", default=False)
parser.add_argument("--membrane", action="store_true", default=False)
parser.add_argument("--hybrid", action="store_true", default=False)


args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

proteinName = pdb_id = args.protein
# chain='A'
# chain='ABC'
chain = args.chain.upper()

pdb = f"{pdb_id}.pdb"

# print(args)
#with open('create_project_commandline_args.txt', 'w') as f:
#    f.write(' '.join(sys.argv))
#    f.write('\n')

# Download the file and rename it to crystal_structure.pdb
#if not os.path.exists(f"crystal_structure.pdb"):
#    pdb_list = [proteinName]
#    downloadPdb(pdb_list)
#    cleanPdb(pdb_list, chain="-1", toFolder="cleaned_pdbs")
#    do(f"cp cleaned_pdbs/{pdb} crystal_structure.pdb")
#
#
#
if chain == "-1":
    chain = getAllChains(pdb)
    print("Chains to simulate: ", chain)
#
## for compute Q
#input_pdb_filename, cleaned_pdb_filename = prepare_pdb("crystal_structure.pdb", chain)
#ensure_atom_order(input_pdb_filename)
# get fasta, pdb, seq file ready
getSeqFromCleanPdb(pdb, chains=chain, writeFastaFile=True)
#do(f"cp crystal_structure.fasta {pdb_id}.fasta")

#if args.extended:
#    do(f"{OPENAWSEM_LOCATION}/helperFunctions/fasta2pdb.py "+proteinName)
##    add_chain_to_pymol_pdb(pdb)  # only work for one chain only now
#else:
#    do(f"cp crystal_structure.pdb {pdb}")

input_pdb_filename, cleaned_pdb_filename = prepare_pdb(pdb, chain)

#ensure_atom_order(input_pdb_filename)

#os.system(f"cp {OPENAWSEM_LOCATION}parameters/burial_gamma.dat .")
#os.system(f"cp {OPENAWSEM_LOCATION}parameters/gamma.dat .")
#os.system(f"cp {OPENAWSEM_LOCATION}parameters/membrane_gamma.dat .")
#os.system(f"cp {OPENAWSEM_LOCATION}parameters/anti_* .")
#os.system(f"cp {OPENAWSEM_LOCATION}parameters/para_* .")
#
##do(f"python {OPENAWSEM_LOCATION}/helperFunctions/Pdb2Gro.py crystal_structure.pdb amh-go.gro")

## ssweight
#do("stride crystal_structure.pdb > ssweight.stride")
#do(f"python {OPENAWSEM_LOCATION}/helperFunctions/stride2ssweight.py > ssweight")

# below used for zim and zimPosition file
#if args.membrane or args.hybrid:
#    do("grep -E 'CB|CA  GLY' crystal_structure-cleaned.pdb > cbs.data")
#    do("""awk '{if($9>15) print "1"; else if($9<-15) print "3"; else print "2"}'  cbs.data  > zimPosition""")
#    do("python3 ~/opt/create_zim.py")
#
#
#if args.frag:
#    do(f"cp crystal_structure.fasta {proteinName}.fasta")
#    do(f"cp {OPENAWSEM_LOCATION}/database/cullpdb_pc80_* .")
#    do(f"python {OPENAWSEM_LOCATION}/helperFunctions/MultCha_prepFrags_index.py \
#    cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 1 9 > logfile" % proteinName)
#    check_and_correct_fragment_memory("frags.mem")
#
#do(f"cp {OPENAWSEM_LOCATION}mm_run.py .")
#do(f"cp {OPENAWSEM_LOCATION}mm_analysis.py .")
