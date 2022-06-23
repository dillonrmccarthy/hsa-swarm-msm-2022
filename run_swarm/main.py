#!/usr/bin/python3
import schrodinger
from schrodinger.structure import StructureReader
from schrodinger.structure import StructureWriter
from schrodinger.application.desmond import cms
from schrodinger.application.desmond.packages import topo
from schrodinger.application.desmond.packages import analysis
from schrodinger.application.desmond.packages import traj
from write_schrodfiles import writecms, writecfg, writemsj, writemsj_norelax
from extrastuff import Bcolors, printdots, gen_seeds
from interpret_cvs import Aall # import for returning top frames and corresponding repelica's
import numpy as np
import os
import sys
import subprocess

#===============================================================================
# the purpose of this script is to:
#   1. Process a generation of 5 replicas (using bash script and VMD)
#   2. analyze the collective variables (using interpret_cvs import)
#       2.1. it analyzes the top 5 of each replica, but returns the top 5 of all replica's
#       2.2. the two lists are both arrays and correspond to eachother. see interpret_cvs.py for more details
#   3. Create a new five new CMS files based off these coordinates after using the "template.cms" file
#   4. Generate a new MSJ and CFG that matches the new CMS, and add them a new replica folder, in a new generation.
#
# So basically this script analyzes the current generation's replica's and outputs a new generation based off the collective variables.
#
# THIS SCRIPT SHOULD BE RUN WITHIN THE CURRENT GENERATION, at the replica level directory.
#
#===============================================================================
if __name__ == "__main__":
    _rcheck = len([fh for fh in sorted(os.listdir()) if "rep_" in fh])
    if _rcheck < 5:
        print("you are not in the correct directory")
        exit()

#-------------------------------------------------------------------------------
#NOTE: defining the core file-paths, and other key information
while True:
    cdir = os.path.abspath(os.getcwd()) #path of where script was submitted from
    pdir = os.path.abspath('..')
    script_dir = '/home/dillon/Scripts/dna_hsa/shotgun_sim/1_shotgun' #directory where ALL scripts are located
    template = os.path.join(pdir,"template.cms") #template cms file path
    bash_submit_main = os.path.join(script_dir,'submit_vmd_main.sh') # file path for the bash submit script
    tcl_script = os.path.join(script_dir,'generate_cvs.tcl') #file path for the tcl script, which is what is passed into VMD in the actual bash script
    vmd_indc = os.path.join(cdir,'.VMD_COMPLETE.INDC')

    #loading indices arrays, for indexing specific subsets of the system
    solu_p = os.path.join(pdir, 'solute_indc.npy')
    sod_p = os.path.join(pdir, 'sodium_indc.npy')
    solv_p = os.path.join(pdir, 'solvent_indc.npy')

    solu_indx = np.load(solu_p)
    sod_indx = np.load(sod_p)
    solv_indx = np.load(solv_p)

#    print(len(solu_indx))
#    print(len(sod_indx))
#    print(len(solv_indx))
#    print(len(solu_indx)+len(sod_indx)+len(solv_indx))
    break

#-------------------------------------------------------------------------------
#read in a trajectory, return the correct xyz coords and box size for a certain frame
def get_new_prop(cms_file, traj_path,frame):
    _throwaway = cms_file
    #msys_modelt, cms_modelt = topo.read_cms(cms_file) #this makes everything INSANELY slow so dont use it
    trj = traj.read_traj(traj_path)
    trj2 = traj.Source(trj).retrieve_frame(frame)
    xyz = trj2.pos()
    box = trj2.box[0,0]
    return xyz, box

#-------------------------------------------------------------------------------
#make new generation and replicas
def make_new_gen(pdir, gen):
    gen_new = "gen_"+str((int(gen)+1))
    _ngp = os.path.join(pdir,gen_new)
    if len([fh for fh in os.listdir(pdir) if gen_new in fh]):
        print("CHECK CODE, YOU MAY BE OVERWRITING A PRIOR GENERATION")
        sanity_chk = input("DO YOU WANT TO CONTINUE? THIS MAY ERASE EVERYTHING!(y/N): ")
        if sanity_chk == 'y' or sanity_chk == 'Y':
            pass
        else:
            exit()
            pass
    else:
        os.mkdir(_ngp)
    return _ngp

#-------------------------------------------------------------------------------
#make the new replica folders within the new generation folder
def make_new_replica_folder(new_gen_path, rep_num, sim_folder):
    rep_new = "rep_"+str(rep_num)
    _nrp = os.path.join(new_gen_path,rep_new)
    _nsp = os.path.join(_nrp,sim_folder)
    if len([fh for fh in os.listdir(new_gen_path) if rep_new in fh]):
        #print("CHECK CODE, YOU MAY BE OVERWRITING A PRIOR REPLICA")
        #sanity_chk = input("DO YOU WANT TO CONTINUE? THIS MAY ERASE EVERYTHING!(y/N): ")
        sanity_chk= 'y'
        if sanity_chk == 'y' or sanity_chk == 'Y':
            pass
        else:
            exit()
            pass
    else:
        os.mkdir(_nrp)
        os.mkdir(_nsp)
    return _nsp

#-------------------------------------------------------------------------------
#take indicies, master xyz array, return new xyz array...
def make_subset_cords(cords): #cords is main xyz, loi is list of indices for whatever
    solute_i = solu_indx-1
    sodium_i = sod_indx-1
    solvent_i = solv_indx-1

    #solute, sodium, solvent = cords[solute_i], cords[sodium_i], cords[solvent_i]
    #num_solute, num_sodium, num_solvent = len(solute_i), len(sodium_i), len(solvent_i)

    xyz_dict = {
        "solute_xyz" : cords[solute_i],
        "solute_num" : len(solute_i),
        "sodium_xyz" : cords[sodium_i],
        "sodium_num" : len(sodium_i),
        "solvent_xyz" : cords[solvent_i],
        "solvent_num" : len(solvent_i),
        "all_num" : (len(solu_indx)+len(sod_indx)+len(solv_indx)) #total number of atoms
    }
    return xyz_dict

#-------------------------------------------------------------------------------
#First, check to see if the VMD analysis has been done, run analysis if it has not
file_exists = os.path.exists(vmd_indc) #does the indicator file exist?
if not file_exists:
    #_vmdskip = 'yes'
    _vmdskip = 'no'
    subprocess.run(["sh",bash_submit_main,tcl_script,cdir,script_dir,_vmdskip]) #run the vmd processing if it hasnt been run
else: printdots(f"{Bcolors.HEADER}Skipping VMD Analysis")
print(f"{Bcolors.ENDC} ") #reset color


#and also check which generation we are on
current_generation = str(cdir[-1:])
gen_check_catch = input(f"Analysing generation {current_generation}. Is this correct? (Y/n): ".format(current_generation))
if not gen_check_catch or gen_check_catch == "y" or gen_check_catch == "Y": pass
else: print(f"{Bcolors.FAIL}CHECK GENERATION{Bcolors.ENDC}"); exit()

#-------------------------------------------------------------------------------
#analyze all the collective variables, return the top 5 frames, using interpret_cvs.py
frames,replicas = Aall()
new_gen_path = make_new_gen(pdir,current_generation) #create path for new generation and makes the folder
seeds = gen_seeds(current_generation)
print("\n")
print(f"{Bcolors.b}{Bcolors.u}Top Frames are:{Bcolors.ur} %s{Bcolors.ENDC}" % (frames))
print(f"{Bcolors.b}{Bcolors.u}Top Replicas are:{Bcolors.ur} %s{Bcolors.ENDC}" % (replicas))
print("\n")

for n, (replica, frame) in enumerate(zip(replicas,frames)):
    print("working on frame",str(frame),"replica",str(replica))
    print(" ")

    seed = seeds[n]
    g_info,r_info,f_info = "g-"+current_generation,"r-"+str(replica),"f-"+str(frame) #this is all the infor for which frame i need to pull!
    base_name = g_info+"_"+r_info+"_"+f_info #this will be the base name of the next iteration of the simulation
    rep_folder = "rep_"+str(replica)
    rep_folder_path = os.path.join(cdir,rep_folder) #gives the file path of the replica that contains the frame we need

    #get all the info we need
    sim_folder = os.listdir(rep_folder_path)[0] #this is one down from the parent replica folder
    sim_folder_path = os.path.join(rep_folder_path,sim_folder) #this is the path of the sim folder main
    traj_folder = [xxx for xxx in os.listdir(sim_folder_path) if "_trj" in xxx][0] #this is the name of the trajectory folder
    traj_path = os.path.join(sim_folder_path, traj_folder) #this is the path of the actual traj file

    #now get the info for those frames
    cord, boxsize = get_new_prop(None, traj_path,frame) #dont need to give the cms file, only trj

    new_replica_path = make_new_replica_folder(new_gen_path, n, base_name)
    cms_fn, msj_fn, cfg_fn  = base_name+".cms",base_name+".msj",base_name+".cfg"
    cms_file_path = os.path.join(new_replica_path,cms_fn)
    msj_file_path = os.path.join(new_replica_path,msj_fn)
    cfg_file_path = os.path.join(new_replica_path,cfg_fn)

    xyz_dict = make_subset_cords(cord)
    xyz_dict['all_xyz'] = cord #add master cords to dict

    writecms(template, xyz_dict, boxsize,cms_file_path) #write cms
    #writemsj(msj_file_path, cfg_fn)
    writemsj_norelax(msj_file_path, cfg_fn)
    writecfg(cfg_file_path, seed)

exit()

#Old Code
"""
object_methods = [method_name for method_name in dir(st) if not callable(getattr(st, method_name))]
print(object_methods)
#print(st.box)

#in_cms = [fh for fh in sorted(os.listdir()) if fh[-4:] == ".cms" and "in" not in fh and "out" not in fh][0]

#printing with the f thing
current_generation = int(cdir[-1:])
#gen_check_catch = input(f"Analysing generation {current_generation:6}. Is this correct?".format(current_generation))
"""
