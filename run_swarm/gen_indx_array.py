import schrodinger
from schrodinger.structure import StructureReader
from schrodinger.structure import StructureWriter
from schrodinger.application.desmond import cms
from schrodinger.application.desmond.packages import topo
from schrodinger.application.desmond.packages import analysis
from schrodinger.application.desmond.packages import traj
from write_schrodfiles import writecms, writecfg, writemsj
from interpret_cvs import Aall # import for returning top frames and corresponding repelica's
import numpy as np
import os
import sys
import subprocess

import pickle as pkl

if __name__ == '__main__':
    print("no")
    #sys.exit(1)

def gen_indx_arrays(cms_file, traj_path,frame):
    sys_cms = cms.Cms(cms_file)
    #define asl for subtypes
    solute_asl = '(( all) AND NOT ((res.ptype "SPC "))) AND NOT ((res.ptype "NA  "))'
    sod_asl = '(res.ptype "NA")'
    solvent_asl = '(res.ptype "SPC ")'
    #y = x.select_atom_comp(asl)
    x = np.array(sys_cms.select_atom(solute_asl))
    y = np.array(sys_cms.select_atom(sod_asl))
    z = np.array(sys_cms.select_atom(solvent_asl))
    np.save("solute_indc",x)
    np.save("sodium_indc",y)
    np.save("solvent_indc",z)

    #msys_modelt, cms_modelt = topo.read_cms(cms_file) #this makes everything INSANELY slow so dont use it
    #trj = traj.read_traj(traj_path)
    ##print(trj[0].pos())
    #tmp = traj.Source(trj).__dict__
    #trj2 = traj.Source(trj).retrieve_frame(frame)
    ##print(trj2.pos())
    ##print(trj2.box)
    #xyz = trj2.pos(z)
    #box = trj2.box[0,0]
    #print(xyz)
    #return xyz, box

#traj_path = '/home/dillon/HDD_C/HSA_Project/simulations/shotgun_simulation/c4_in_main/gen_0/rep_0/c4in-v2_f-10416_a/c4in-v2_f-10416_a_trj'
#frame = 0
#cms_file = '/home/dillon/HDD_C/HSA_Project/simulations/shotgun_simulation/c4_in_main/gen_0/rep_0/c4in-v2_f-10416_a/c4in-v2_f-10416_a.cms'
#cms_file = '/home/dillon/HDD_C/HSA_Project/simulations/shotgun_simulation/c4_in_main/template2.cms'
#gen_indx_arrays(cms_file,None,None)

s1 = np.load('sodium_indc.npy')
s2 = np.load('/home/dillon/HDD_C/HSA_Project/simulations/shotgun_simulation/c4_in_main/gen_0/rep_0/c4in-v2_f-10416_a/sodium_indc.npy')
t1 = np.load('solute_indc.npy')
t2 = np.load('/home/dillon/HDD_C/HSA_Project/simulations/shotgun_simulation/c4_in_main/gen_0/rep_0/c4in-v2_f-10416_a/solute_indc.npy')
u1 = np.load('solvent_indc.npy')
u2 = np.load('/home/dillon/HDD_C/HSA_Project/simulations/shotgun_simulation/c4_in_main/gen_0/rep_0/c4in-v2_f-10416_a/solvent_indc.npy')
#print(np.array_equal(s1,s2))
#print(np.array_equal(t1,t2))
#print(np.array_equal(u1,u2))

#print(t1)
