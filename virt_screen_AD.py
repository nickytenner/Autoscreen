import argparse
import numpy as np
import os
import tempfile
import pyprind


def main(parser):
    #parser.add_argument("positional1",help="Add a positional argument like this if needed")
    parser.Named = parser.add_argument_group("required arguments")
    parser.Named.add_argument("--irec", help="Receptor pdb file", type=str, action="store", required=True)
    parser.Named.add_argument("--dlig", help="Directory containing ligand sdf files", type=str, action="store", required=True)
    parser.add_argument("--sdflist", default=None, help="List of screened sdfs obtained using MLMM", type=str, action="store")
    parser.add_argument("--prob", default=1.0, help="Enter the fraction of ligands to be randomly selected from the sdf folder", type=float, action="store")
    #parser.add_argument("--odpf",default="dock.dpf",help="Docking parameter file to be written",type=str,action="store")
    # parser.add_argument("--grid",default = [100,100,100],\
    #                    help="Enter the grid points as x y z",type=int,nargs="+",action="store")
    # parser.Named.add_argument("--center",default = [100,100,100],\
    #                          help="Enter the center of the grid box as x y z",type=float,nargs="+",action="store",required=True)
    # parser.add_argument("--spacing",default=0.2,type=float,action="store")
    args = parser.parse_args()
    #print (args.ipdbqt,args.ogpf,args.grid)
    receptor = args.irec
    ligdir = args.dlig
    fsdflist = args.sdflist
    prob = args.prob
    sdflist = []
    if fsdflist:
        with open(fsdflist) as fp:
            fp.readline()
            for line in fp:
                sdflist.append(line.split()[1])
        print ("Requested {} sdf files via {}".format(len(sdflist), fsdflist))
    liglist = []
    # print(receptor,ligdir)
    print ("Starting virtual screening using AutoDock 4.1")
    print ("Receptor:" + receptor)
    print ("Ligand directory:" + ligdir)
    for file in os.listdir(ligdir):
        if file.endswith(".sdf"):
            if fsdflist:
                if file in sdflist:
                    liglist.append(file)
            else:
                liglist.append(file)
    print (str(len(liglist)) + " ligands found")
    assert prob <= 1.0 and prob >= 0.0, "ERROR: Fraction must be between 0 and 1 !!"
    # if nmax == 0:
    #    prob = 1.0
    # else:
    #    prob = float(nmax)/(float(len(liglist)))
    if fsdflist:
        if len(sdflist) != len(liglist):
            print ("WARNING!! Not all specified sdfs ({}) were found ({}) in the ligand folder.".format(len(sdflist), len(liglist)))
            for elem in sdflist:
                if elem not in liglist:
                    print ("Ligand {} not found in the folder {}... Continuing".format(elem, ligdir))
    #print (liglist)
    # quit()
    print ("tail logfile to track the run")
    try:
        os.system('rm logfile')
    except:
        pass
    bar = pyprind.ProgBar(len(liglist), bar_char='#')
    bar.update()
    docking_stat = {}
    os.system('cp ' + receptor + ' receptor.pdb')
    os.system('sh ./scripts/prep_recept.job')
    for lig in liglist:
        if np.random.rand() < prob:
            rundir = tempfile.TemporaryDirectory(dir="./")
            # print(rundir.name)
            os.system('cp receptor.pdbqt ' + str(rundir.name) + '/')
            os.system('cp ' + os.path.join(ligdir, lig) + ' ' + str(rundir.name) + '/ligand.sdf')
            os.system('cp ' + "scripts/AD4.1_bound.dat" + ' ' + str(rundir.name))
            os.system('cp ' + "scripts/dock.job" + ' ' + str(rundir.name))
            #os.system('ls '+str(rundir.name))
            os.chdir(str(rundir.name))
            # os.system('pwd')
            # os.system('ls')
            os.system('sh ./dock.job')
            dlg_file = receptor[:-4] + '_' + lig[:-3] + 'dlg'
            xyz_file = receptor[:-4] + '_' + lig[:-3] + 'xyz'
            os.system('cp dock.dlg ../' + dlg_file)
            os.system('cp ligand.xyz ../' + xyz_file)
            os.chdir('../')
            rundir.cleanup()
            with open(dlg_file) as fp:
                for line in fp:
                    if "Cluster Rank = 1" in line:
                        for i in range(5):
                            en_info = fp.readline()
                        docking_stat[lig] = float(en_info.split()[7])
                        break
        bar.update()
    with open("logfile", "a") as fp:
        print ("Docking results")
        fp.write("Docking results\nRank; Ligand ; Score\n")
        lig_cand = sorted(docking_stat, key=docking_stat.get)
        for i, lig in enumerate(lig_cand):
            print (lig, " : ", docking_stat[lig], "kcal/mol")
            fp.write("%d,  %s , %6.2f kcal/mol\n" % (i, lig, docking_stat[lig]))
        print("NORMAL TERMINATION")
        fp.write("NORMAL TERMINATION")


if __name__ == "__main__":
    main(argparse.ArgumentParser())
