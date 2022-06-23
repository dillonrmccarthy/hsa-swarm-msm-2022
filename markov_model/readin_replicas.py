#!/usr/bin/python3

import os
import numpy as np

def gen_paths(cv_fname,**kwargs): #creates a list of all the data files, in order of gen/replica (see readme.md)
    #all_gens = sorted([x for x in os.listdir() if x[0:4] == "gen_"],key=lambda y:int(y[5:]))
    #all_gens =[y for y in os.listdir() if y[0:4] == "gen_"]
    #sorted_all_gens = sorted(all_gens, key=lambda x: int(x[4:]))
    #sorted_gens = sorted([y for y in os.listdir() if y[0:4] == "gen_"], key=lambda x: int(x[4:]))
    #print(os.path.abspath(sorted_gens[0])) #returns the generation abs path
    #full_gen_paths = [os.path.abspath(gn) for gn in sorted_gens]

    #sweaty one liner
    if 'nreps' in kwargs:
        _nreps = int(kwargs['nreps'])
    else:
        _nreps = False

    full_gen_paths = [os.path.abspath(gn) for gn in sorted([y for y in os.listdir() if y[0:4] == "gen_"], key=lambda x: int(x[4:]))]
    #print(len(full_gen_paths))

    """
    ok now making a really interesting array
    the shape is (num of generations, number of replicas, [information])
    so...
    [[[gen 0, rep 0, o_gen, o_rep, o_frame, path-to-dat ... (basically reading in the first data file)

    for now, the information will JUST be the abs file path. this gets fed to make LOA
    """
    all_data_file_paths_1 = []
    all_data_file_paths_2 = []
    for gen in full_gen_paths:
        td = 0
        _gnme= os.path.basename(gen)
        for root, dirs, files in os.walk(gen,topdown=True):
            dirs[:] = [d for d in dirs if d[0:4] =="rep_"] #filter out to only directories with ht_rep in them BE CAREFUL

            #sorted_rep = [y for y in [os.path.join(root,rep_name) for rep_name in sorted(dirs,key=lambda x:int(x[4:]))] if len(y)>1] # this is a dumb way tbh
            if len(dirs)>0:
                sorted_rep = [os.path.join(root,rep_name) for rep_name in sorted(dirs,key=lambda x:int(x[4:]))] #retuns a list of all the rep's in order g0(r0,r1...) g1(r0,r1..)
                #will have to os.walk this list to get the actual dat files.
                for rep_folder in sorted_rep:
                    for root, dirs, files in os.walk(rep_folder):
                        if cv_fname in files:
                            datfile = [dfile for dfile in files if dfile == cv_fname]
                            #print(root,datfile[0])
                            all_data_file_paths_1.append(os.path.join(root,datfile[0]))

                        #doing it twice just to be safe
                        for _file in files:
                            tmp_list = []
                            printer = 0
                            if _file == cv_fname: #hey look the input is finally back!
                                all_data_file_paths_2.append(os.path.join(root,_file))

            else:
                continue


    #sanity check
    if all_data_file_paths_1 != all_data_file_paths_2:
        print("uh oh.....")
        exit()

    #for h in range(len(all_data_file_paths_1)):
    #    print(all_data_file_paths_1[h].split("/")[-4:])
    #    print(all_data_file_paths_2[h].split("/")[-4:])
    #    print("\n\n")
    return all_data_file_paths_1




#NOTE: this is currently structured such that each simulation folder has only one file with data points (i.e, cannot take multiple files of data for a single trajectory. If this is required, add a second function which concatenates data.


#this takes cvp, which is the list of all cv data file paths
class make_loa: #takes a list of all file paths for all cv.dat files, makes list of arrays
#modify later?
    def __init__(self,cvp,**kwargs):
        self.cvd = [] #cvd=cv data master list

        for datf in cvp:
            numpy_array_made = False
            """for each file, we read in a line of [a, b, c ..., n]
            what we want as the output, is an array that looks like this:
                array( [[a1,b1,c1...n1]
                       [a2,b2,c2...n2]
                       [a2,b2,c2...n2]] ) where the shape is (number of frames, cvs per frame)"""

            with open(datf,'r') as fh:
                for line in fh.readlines():
                    if "#" not in line: #read in the very first line
                        tmp = list(map(float,line.replace(" ","").replace("\n","").split(',')))
                    if numpy_array_made == False:
                        try:
                            tmp
                        except UnboundLocalError:
                            continue #if does not exist, continue
                        else: #if it does exist
                            traj_ary = np.array(tmp) #make your initial array [3,] array
                            traj_ary = np.expand_dims(traj_ary,axis=0) #expand so you have a [1,3]
                            numpy_array_made = True #make sure you cant come back here
                            continue #continue to the next line because you dont want to hit the next if statement

                    if numpy_array_made == True:
                        tmp_ary = np.array(tmp) # gives array of shape (len-of-tmp,)
                        #we need it to be an array of (1,len) so that we can concatenate it with any shape (n,len) array along axis=0
                        tmp_ary = np.expand_dims(tmp_ary,axis=0)
                        traj_ary = np.concatenate([traj_ary,tmp_ary],axis=0)


                        #finally:
                            #print("---")
                fh.close()
            self.cvd.append(traj_ary)

    def make_1d(self,n):# n is the dimension you want to slice along
        loa = self.cvd
        for i,ind_traj in enumerate(loa):
            if i == 0:
                m_array = ind_traj[:,n].T #just grabbing the distance measuremen
            else:
                single = ind_traj[:,n].T #just grabbing the distance measuremen
                m_array = np.concatenate([m_array,single],axis=0)

        return m_array


    def strip_dim(self,n,m): #this strips demensions off the arrays in the big list
        m +=1
        new_array = []
        loa = self.cvd
        for orig in loa:
            tmp = orig[:,n:m]
            new_array.append(tmp)

        return new_array



        self.definitions = {
        'cvd':'stands for Collective Variable Data. A list of arrays, each array is 1 traj of collective variabless'#,
        }


if __name__ == '__main__':
    exit()
