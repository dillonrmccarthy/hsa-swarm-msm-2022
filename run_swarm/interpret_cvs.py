#!/usr/bin/python3

#This interprets the CV's of all 5 replicas, and returns the top 5 frames from all replicas as two numpy arrays, where array 1 = the top 5 frames [1,2,...5] and array 2 = the replica's those frames came from [1,1,...4]

'''
example return data:
[941,332,615,346,222] <--- frames
[  3,  2,  2,  2,  4] <--- replicas

frame 941 of replica 3, frame 332 of replica 2, etc...
'''

import numpy as np
import sys, os
from textwrap import wrap
from os import listdir, getcwd
from glob import glob

#-------------------------------------------------------------------------------
class Asr:
    def __init__(self,**kwargs):
        if 'replica' in kwargs.keys():
            self.replica = int(kwargs['replica'])
            #print("Replica:",self.replica)
        else:
            print("need replica")
            exit()

    def analyze(self,filename,**kwargs):
        fh = open(filename,'r')
        m_list = []
        if 'cutoff' in kwargs.keys():
            cutoff = float(kwargs['cutoff'])
        for line in fh.readlines():
            if "#" in line:
                continue
            else:
                tmp = line.replace(" ","").replace("\n","").split(',')
                #loser = lambda : float(kwargs['cutoff']) if 'cutoff' in kwargs.keys() else False
                #if loser():
                if 'cutoff' in kwargs.keys() and float(tmp[2]) >= cutoff:
                    continue
                else:
                    m_list.append([int(tmp[0]), float(tmp[1]), float(tmp[2])])
        tmp2 = np.asarray(m_list)
        self.sorted_frames = tmp2[tmp2[:,1].argsort()]
        self.top5 = self.sorted_frames[-1:-6:-1] ;     #reverse sorted frames
        self.repID = self.replica * np.ones_like(self.top5[:,0])
        while False:
            #self.top5_b = self.sorted_frames[-1:-6:-1][:,0:2] #this gives the reverse sorted, but only first two values
            #self.top5_fr = self.top5[:,0]
            continue





#---------------------------------------------------------------------------------------------

#r = [0,1,2,3,4]
#data = Asr(replica=r[0])
#data.analyze('/home/dillon/Scripts/dna_hsa/shotgun_sim/test.dat',cutoff=30)
#print(data.top5)
#print(data.repID)

#initialize all replica folder paths, and all replica's ID's
def Aall():
    all_cont = listdir(getcwd())
    rep_folders = [j for i,j in enumerate(all_cont) if "rep_" in j]
    for k,m in enumerate(rep_folders):
        replica_path = os.getcwd()+'/'+m
        path = replica_path
        #get subfolder
        replica = int(m.split('_')[-1])
        for root, dirs, files in os.walk(path):
            for file in files:
                if file == 'cv_results.dat':
                #if file == 'test.dat':
                    cv_data = os.path.join(root,file)
                    #print("Working on data file:", cv_data)
                    #print("in simulation folder:",root)
                    data = Asr(replica=replica)
                    #data.analyze(cv_data,cutoff=30)
                    data.analyze(cv_data)
                    #print("TOP 5 FRAMES")
                    #print(data.top5)
                    #print(data.repID)
                    #print(" ")
                    #print(" ")
                    try:
                        cv_results
                    except NameError:
                        cv_results = data.top5
                        replica_ids = data.repID
                    else:
                        cv_results = np.concatenate((cv_results,data.top5),axis=0)
                        replica_ids = np.concatenate((replica_ids,data.repID),axis=0)

    #print("=============================================================")
    #return cv_results, replica_ids
    #print(cv_results)
    #print(replica_ids)

    ind = np.argsort(cv_results[:,1])
    #print(cv_results[ind])
    #print(cv_results[ind[-1:-6:-1]]) #gets the last 6 indicies
    #print(replica_ids[ind[-1:-6:-1]]) #gets the last 6 indicies
    top5_all = cv_results[ind[-1:-6:-1]]
    top5_frames = cv_results[ind[-1:-6:-1],0]
    top5_rep = replica_ids[ind[-1:-6:-1]]
    #vfunc = np.vectorize(int)
    #tmp_t5f = vfunc(top5_frames)
    #tmp_t5r = vfunc(top5_rep)
    tmp_t5f = top5_frames.astype(int) #better way to cast to int
    tmp_t5r = top5_rep.astype(int)



    #print(tmp_t5f,tmp_t5r)
    #return top5_frames,top5_rep
    return tmp_t5f,tmp_t5r #mapped int to all the values so that they are ints not floats

#frames, replicas = Aall()

#print(frames)
#print(replicas)
