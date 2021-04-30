#pdbqt = 'Chlorhexdine.pdbqt'

from collections import Counter
import argparse
import numpy as np

def main(parser):
    #parser.add_argument("positional1",help="Add a positional argument like this if needed")
    parser.Named = parser.add_argument_group("required arguments")
    parser.Named.add_argument("--irec",help="Receptor pdbqt file",type=str,action="store",required=True)
    parser.Named.add_argument("--ilig",help="Ligand pdbqt file",type=str,action="store",required=True)
    parser.add_argument("--ogpf",default="grid.gpf",help="Grid parameter file to be written",type=str,action="store")
    parser.add_argument("--odpf",default="dock.dpf",help="Docking parameter file to be written",type=str,action="store")
    parser.add_argument("--grid",default = [100,100,100],\
                        help="Enter the grid points as x y z",type=int,nargs="+",action="store")
    parser.Named.add_argument("--center",default = [100,100,100],\
                              help="Enter the center of the grid box as x y z",type=float,nargs="+",action="store",required=True)
    parser.add_argument("--spacing",default=0.2,type=float,action="store")
    args = parser.parse_args()
    #print (args.ipdbqt,args.ogpf,args.grid)
    receptor = args.irec
    ligand = args.ilig
    gpf = args.ogpf
    grid = args.grid
    spacing = args.spacing
    center = args.center
    dpf = args.odpf
    atom_types = []
    mol_center = []
    dock_settings = {
        "autodock_parameter_version":["4.2","# used by autodock to validate parameter set"],
        "parameter_file": ["AD4.1_bound.dat","# parameter library filename"],
        "intelec":[" ","# calculate internal electrostatics"],
        "seed":["pid time","# seeds for random generator"],
        "tran0": ["random","# initial coordinates/A or random"],
        "quaternion0": ["random","# initial orientation"],
        "dihe0": ["random","# initial dihedrals (relative) or random"],
        "torsdof": ["9","# torsional degrees of freedom"],
        "rmstol": ["2.0","# cluster_tolerance/A"],
        "extnrg": ["1000.0","# external grid energy"],
        "e0max":["0.0 10000","# max initial energy; max number of retries"],
        "ga_pop_size":["150","# number of individuals in population"],
        "ga_num_evals":["2500000","# maximum number of energy evaluations"],
        "ga_num_generations":["27000","# maximum number of generations"],
        "ga_elitism":["1","# number of top individuals to survive to next generation"],
        "ga_mutation_rate":["0.02","# rate of gene mutation"],
        "ga_crossover_rate":["0.8","# rate of crossover"],
        "ga_window_size":["10","# "],
        "ga_cauchy_alpha":["0.0","# Alpha parameter of Cauchy distribution"],
        "ga_cauchy_beta":["1.0","# Beta parameter Cauchy distribution"],
        "set_ga":[" ","# set the above parameters for GA or LGA"],
        "sw_max_its":["300","# iterations of Solis & Wets local search"],
        "sw_max_succ":["4","# consecutive successes before changing rho"],
        "sw_max_fail":["4","# consecutive failures before changing rho"],
        "sw_rho":["1.0","# size of local search space to sample"],
        "sw_lb_rho":["0.01","# lower bound on rho"],
        "ls_search_freq":["0.06","# probability of performing local search on individual"],
        "set_psw1":[" ","# set the above pseudo-Solis & Wets parameters"],
        "unbound_model":["bound","# state of unbound ligand"],
        "ga_run":["5","# do this many hybrid GA-LS runs"],
        "analysis":[" ","# perform a ranked cluster analysis"]
    }

    for pdbqt in (receptor,ligand):
        print ("Reading ",pdbqt)
        with open (pdbqt) as fp:
            tot_charge = 0.0
            atoms = []
            xyz = []
            for line in fp:
                if line.strip():
                    data = line.split()
                    if data[0] == "ATOM":
                        #print(data)
                        xyz.append([float(data[-7]),float(data[-6]),float(data[-5])])
                        atoms.append(data[-1])
                        tot_charge += float(data[-2])
        mol_center.append(np.mean(np.array(xyz),axis=0))
        atoms_info = Counter(atoms)
        for key,val in atoms_info.items():
            print(key,":",val)
        print("Atoms types(",len(list(atoms_info.keys())),") and numbers(",sum(list(atoms_info.values())),"):")
        print("Total charge in the molecule ",pdbqt,":",tot_charge)
        atom_types.append(list(atoms_info.keys()))
    #atom_types[1].append("Z")
    #print(atom_types)
    #atom_types_rec = list(atom_types.keys())
    with open(gpf,"w") as fp, open(dpf,"w") as fp2:
        fp.write("npts %3d %3d %3d %s # num.grid points in xyz\n" %(grid[0],grid[1],grid[2]," "*21))
        fp.write("gridfld %s %s # grid data file\n"%(receptor[:-5]+"maps.fld"," "*8))
        fp.write("spacing %f %s # spacing(A)\n"%(spacing," "*21))
        fp.write("receptor_types %s %s # receptor atom types\n"%(" ".join(atom_types[0])," "*6))
        fp.write("ligand_types %s %s # ligand atom types\n"%(" ".join(atom_types[1])," "*10))
        #fp2.write("ligand_types %s %s # ligand atom types\n"%(" ".join(atom_types[1][:-1])," "*10))
        fp2.write("ligand_types %s %s # ligand atom types\n"%(" ".join(atom_types[1])," "*10))
        fp.write("receptor %s %s # macromolecule\n"%(receptor," "*10))
        fp.write("gridcenter %4.1f %4.1f %4.1f %s # xyz-coordinates or auto\n"%(center[0],center[1],center[2]," "*11))
        fp.write("smooth 0.5 %s # store minimum energy w/in rad (A)\n"%(" "*27))
        fp2.write("fld %s %s # grid data file\n"%(receptor[:-5]+"maps.fld"," "*12))
        for elem in atom_types[1]:
            fp.write("map %s %s # atom-specific affinity map\n"%(receptor[:-5]+elem+".map", " "*(16-len(elem))))
            if elem != "Z":
                fp2.write("map %s %s # atom-specific affinity map\n"%(receptor[:-5]+elem+".map", " "*(16-len(elem))))
        fp.write("elecmap %s %s # electrostatic potential map\n"%(receptor[:-5]+"e.map"," "*11))
        fp.write("dsolvmap %s %s # dsolvation potential map\n"%(receptor[:-5]+"d.map"," "*10))
        fp2.write("elecmap %s %s # electrostatic potential map\n"%(receptor[:-5]+"e.map"," "*11))
        fp2.write("dsolvmap %s %s # dsolvation potential map\n"%(receptor[:-5]+"d.map"," "*10))
        fp.write("dielectric -0.1465 %s # <0,AD4 distance-dep.diel;>0, constant\n"%(" "*19))
        fp.write("covalentmap 5.0 1000.0 %6.4f %6.4f %6.4f\n"%(center[0],center[1],center[2]))
        fp2.write("move %s %s # small molecule\n"%(ligand," "*14))
        fp2.write("about %6.4f %6.4f %6.4f %s # small molecule center\n"%(mol_center[1][0],mol_center[1][1],mol_center[1][2]," "*9))
        for key,val in dock_settings.items():
            fp2.write("%s %s %s %s\n"%(key,val[0]," "*(36-len(key)-len(val[0])),val[1]))
        #with open(dpf,"w") as fp:
    #    fp.write("")
    #print(dock_settings)
    
if __name__ == "__main__":
    main(argparse.ArgumentParser())
    
