#pdbqt = 'Chlorhexdine.pdbqt'

import sys,getopt
from collections import Counter

def main(argv):
    pdbqt = ''
    gpf = ''
    try:
        opts,args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print ('pdbqt_check.py -i <file.pdbqt>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('pdbqt_check.py -i <file.pdbqt>')
            sys.exit()
        elif opt in ("-i","--ifile"):
            pdbqt = arg
        elif opt in ("-g","--ifile"):
            gpf = arg
            
    with open (pdbqt) as fp:
        tot_charge = 0.0
        atoms = []
        print (gpf)
        for line in fp:
            if line.strip():
                data = line.split()
                if data[0] == "ATOM":
                    print(data)
                    atoms.append(data[-1])
                    tot_charge += float(data[-2])
    atom_types = Counter(atoms)
    for key,val in atom_types.items():
        print(key,":",val)
    print("Atoms types(",len(list(atom_types.keys())),") and numbers(",sum(list(atom_types.values())),"):")
    print("Total charge in the molecule ",pdbqt,":",tot_charge)

if __name__ == "__main__":
    main(sys.argv[1:])
    
