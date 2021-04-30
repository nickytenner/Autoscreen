#mol2 = 'Chlorhexdine.mol2'

import sys,getopt
from collections import Counter

def main(argv):
    mol2 = ''
    try:
        opts,args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print ('mol2_check.py -i <file.mol2>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('mol2_check.py -i <file.mol2>')
            sys.exit()
        elif opt in ("-i","--ifile"):
            mol2 = arg
            
    with open (mol2) as fp:
        tot_charge = 0.0
        atoms = []
        for line in fp:
            if line.strip():
                data = line.split()
                if len(data) > 7:
                    print(data)
                    atoms.append(data[1])
                    tot_charge += float(data[8])
    atom_types = Counter(atoms)
    for key,val in atom_types.items():
        print(key,":",val)
    print("Atoms types(",len(list(atom_types.keys())),") and numbers(",sum(list(atom_types.values())),"):")
    print("Total charge in the molecule ",mol2,":",tot_charge)
    
    

if __name__ == "__main__":
    main(sys.argv[1:])
    
