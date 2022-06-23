#!/usr/bin/python3

'''
    Script by:      Jacob R Remington
                    Dillon R McCarthy
                    Jianing Li(*)
                    Severin T. Schneebeli(**)

    The results of this script are reproduceable given the datafiles in <datafiles>,
which are the results of the CV analysis for all simulations. the MSM is saved as
<bay_hmsm_lag100_samples500.h5>, cluster centers as <cluster_centers.npy>, and mean
first passage time (MFPT) as <mfpts_lag100_samples500.npy>.

(*) Corresponding author
(**) Corresponding author (secondary)
'''

import numpy as np
import sys
import os
import subprocess
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, MaxNLocator, PercentFormatter
import matplotlib.colors as mcolors
from matplotlib import colors
from matplotlib.colors import Normalize
from matplotlib import cm
import pyemma as pym
from pyemma.coordinates import cluster_kmeans as kmeans, cluster_regspace as regspace, assign_to_centers as as_t_c
import tol_colors as tc
from readin_replicas import *
import scipy.stats as st
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
font = {'family' : 'Arial',
            'size'   : 12}

matplotlib.rc('font', **font)
#config=pym.util._config.Config
##print(config.__dict__)
#config.pyemma_njobs=1
#config.omp_num_threads=1
#exit()
#-------------------------------------------------------------------------------#
#step 1: read in data
dat_array = gen_paths('cv_results.dat',nreps=5) #all the dat files
dat = make_loa(dat_array)#,dims="all") #read it in to an object

def load_stuff(datf):
    dats=[]
    with open(datf,'r') as fh:
        for line in fh:
            if "#" not in line: #read in the very first line
                stuff=line.strip().split(",")
                dats.append(np.array([float(x) for x in stuff[1:3]]))
    return np.stack(dats,axis=0)


#### run mode params
out_directory='<path to filder containing the bay_hmsm_lag100_samples500.h5,cluster_centers.npy' #this is where the MSM out folder is


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#  IF Clustering / MSM / Mean First passage Time is already done, set to 0 and this will reload everything
#
recluster=0#to recluster or not to recluster is the question
recalc_MSM=0#to restimate the bayesian MSM
recalc_mfpt=0

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#this strips the arrays to [(n,2),(n,2)...]
r_loa = dat.strip_dim(1,2) #if you want it to the end you have to do it jenky sorry
r_loa=[x[5:] for x in r_loa]#skip first 5 frames
#
indp_trajs_in=[]
indp_trajs_out=[]
for i in range(1,4):
    stuf=load_stuff("/home/jacob/More_Data/Dillon/shotgun_sim-msm_current/master_batfiles/c4_hsa_out/v"+str(i)+"/cv_results.dat")
    indp_trajs_out.append(stuf)
    stuf=load_stuff("/home/jacob/More_Data/Dillon/shotgun_sim-msm_current/master_batfiles/c4_hsa_in/v"+str(i)+"/cv_results.dat")
    indp_trajs_in.append(stuf)
for i in range(10):
    indp_trajs_in.append(r_loa[i])
r_loa=r_loa[10:]#skip first 2 iterations
r_loa_cat=np.concatenate(r_loa,axis=0)

indp_trajs_in=np.concatenate(indp_trajs_in,axis=0)
indp_trajs_out=np.concatenate(indp_trajs_out,axis=0)

all_data=np.concatenate([indp_trajs_in,r_loa_cat,indp_trajs_out],axis=0)

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#compute relative free energy of independent trajs
kde_in = st.gaussian_kde(indp_trajs_in.T)#
kde_out = st.gaussian_kde(indp_trajs_out.T)
nlevl=100
ngrid=100
minn,maxx=np.min(all_data,axis=0),np.max(all_data,axis=0)
X,Y=np.meshgrid(np.linspace(minn[0],maxx[0],ngrid),np.linspace(minn[1],maxx[1],ngrid))
XY=np.stack([X.flatten(),Y.flatten()],axis=1).T
#compute probability at each grid point
kde_Prob_prediction_in=np.reshape(kde_in(XY),(ngrid,ngrid))
kde_Prob_prediction_out=np.reshape(kde_out(XY),(ngrid,ngrid))
#convert to free energy
#ensure sum ->1
kde_Prob_prediction_in=kde_Prob_prediction_in/np.sum(kde_Prob_prediction_in)
kde_Prob_prediction_out=kde_Prob_prediction_out/np.sum(kde_Prob_prediction_out)
kT=1.987*(300)/1000#kcal/molK assumes
free_energy_in=-np.log(kde_Prob_prediction_in)*kT
free_energy_out=-np.log(kde_Prob_prediction_out)*kT
free_energy_in-=np.min(free_energy_in)
free_energy_out-=np.min(free_energy_out)

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#KMEANS TIME!
if recluster:
    clust_obj_1 = kmeans(r_loa,k=50,max_iter=1000,n_jobs=1)
    #clust_obj_1 = regspace(r_loa,dmin=,max_iter=1000,n_jobs=1)# this is for looking at the regspace
    clust_cent_1 = clust_obj_1.clustercenters
    dtrajs_1=clust_obj_1.dtrajs
    np.save(out_directory+"cluster_centers.npy",clust_cent_1)

else: #this reloads the cluster centers!
    clust_cent_1=np.load(out_directory+"cluster_centers.npy")
    dtrajs_1=ass_t_c(r_loa,clust_cent_1)


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#Plot the cluster centers
norm=matplotlib.colors.LogNorm()
figure=plt.figure(figsize=(3.5,3.5))
plt.hist2d(np.concatenate(r_loa,axis=0)[:,0],np.concatenate(r_loa,axis=0)[:,1],norm=norm,bins=20,cmap='nipy_spectral')
plt.tick_params(axis='both',direction='in')
plt.xticks([50,60,70,80])
plt.yticks([0.,5.,10,15,20])
plt.xlabel("Out of Plane Distance / "+r"$\AA$")
plt.ylabel("In Plane Distance / "+r"$\AA$")
cbar2=plt.colorbar()
cbar2.ax.get_yaxis().labelpad = 15
cbar2.set_label('Counts', rotation=270)
plt.scatter(clust_cent_1[:,0],clust_cent_1[:,1],color='k')
plt.tight_layout()
plt.savefig(out_directory+"Figure_clusters.png",transparent=True,dpi=300)
plt.show()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#Plot the implied timescales
tica_imp_obj_1=pym.msm.its(dtrajs_1,lags=np.arange(1,300,5),reversible=True)#,errors='bayes',nsamples=50)

pym.plots.plot_implied_timescales(tica_imp_obj_1,ylog=True,nits=10,units='ns',dt=0.01)#units='steps')
plt.savefig(out_directory+"Figure_impts.png",transparent=True,dpi=300)
plt.show()


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# MSM Parameters
lag = 100 #chosen from plot of implied timescales
timestep = 10 # 10ps recording inverval
convert = 1000/(lag*timestep)
convert2=1000/timestep
nsampl=500
print(convert)
# x transitions/lag_time * 1 lag/30 steps * 1 step/20 ps * 1000 ps/1 ns

#Estimate the MSM
#Normal MSM
if recalc_MSM:
    msm_1 = pym.msm.bayesian_markov_model(dtrajs_1, lag, reversible=True,nsamples=nsampl) #this is the msm itself
    msm_1.save(out_directory+'bay_hmsm_lag'+str(lag)+'_samples'+str(nsampl)+'.h5',overwrite=True)
else:
    msm_1=pym.load(out_directory+'bay_hmsm_lag'+str(lag)+'_samples'+str(nsampl)+'.h5') #loads in the MSM if already done

#plot slowest eigenvector
def plot_veronin_cmapped(axed,vec_to_color):
    #colorrs=(vec_to_color)
    normalize = matplotlib.colors.Normalize(vmin=np.min(vec_to_color), vmax=np.max(vec_to_color))
    cmap=plt.get_cmap("bwr")
    scalarMap = matplotlib.cm.ScalarMappable(norm=normalize, cmap=cmap)
    points=clust_cent_1[msm_1.active_set]
    big_n=100.
    big_points=np.array([[-big_n,-big_n],[-big_n,big_n],[big_n,-big_n],[big_n,big_n]])
    points=np.concatenate([points,big_points],axis=0)
    vor = Voronoi(points)

    voronoi_plot_2d(vor,ax=axed,show_points=False,show_vertices=False,point_color='k')
    for region in vor.regions:
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            if len(polygon)>0:
                polygon_shapely = Polygon(polygon)
                #determine which state this is
                for activei,statei in enumerate(msm_1.active_set):
                    point = Point(clust_cent_1[statei,0],clust_cent_1[statei,1])
                    if polygon_shapely.contains(point):
                        ind=activei
                        break
                axed.fill(*zip(*polygon),color=scalarMap.to_rgba(vec_to_color[ind]))#,cmap=the_cmap)#plt.get_cmap("bwr")(vec_to_color[ind]),vmin=np.min(vec_to_color),vmax=np.max(vec_to_color)
    axed.set_xlim([np.min(r_loa_cat[:,0]),np.max(r_loa_cat[:,0])])
    axed.set_ylim([np.min(r_loa_cat[:,1]),np.max(r_loa_cat[:,1])])
    axed.tick_params(axis='both',direction='in')
    axed.set_xlabel("Out of Plane Distance / "+r"$\AA$")
    axed.set_ylabel("In Plane Distance / "+r"$\AA$")
    axed.scatter(clust_cent_1[msm_1.active_set,0],clust_cent_1[msm_1.active_set,1],c='k',alpha=1,s=3, zorder=10)

fig=plt.figure(figsize=(3.5,3.5))
axe=fig.add_subplot()
plot_veronin_cmapped(axe,msm_1.eigenvectors_right()[:,i])
plt.tight_layout()
plt.savefig(out_directory+"Figure_slowest_evec.png",transparent=True,dpi=300)
plt.show()



#may need to account for the active_set if pyemma removed some states
kde = st.gaussian_kde(clust_cent_1[msm_1.active_set].T,weights=msm_1.pi)#,bw_method=0.5)#
#minn,maxx=np.min(r_loa_cat,axis=0),np.max(r_loa_cat,axis=0)
#X,Y=np.meshgrid(np.linspace(minn[0],maxx[0],ngrid),np.linspace(minn[1],maxx[1],ngrid))
#XY=np.stack([X.flatten(),Y.flatten()],axis=1).T
#compute probability at each grid point
kde_Prob_prediction=np.reshape(kde(XY),(ngrid,ngrid))
#convert to free energy
#ensure sum ->1
kde_Prob_prediction=kde_Prob_prediction/np.sum(kde_Prob_prediction)

free_energy_msm=-np.log(kde_Prob_prediction)*kT
free_energy_msm-=np.min(free_energy_msm)


fig,axes=plt.subplots(ncols=3,figsize=(6.5,3))

energy_cut=10.
free_energy_in[free_energy_in>=energy_cut]=energy_cut
free_energy_out[free_energy_out>=energy_cut]=energy_cut
free_energy_msm[free_energy_msm>=energy_cut]=energy_cut


pltt=axes[0].contourf(X,Y,free_energy_in,levels=nlevl,cmap='nipy_spectral')
cbar0=plt.colorbar(pltt,ax=axes[0])
cbar0.ax.get_yaxis().labelpad = 15
axes[0].set_xlabel("Out of Plane Distance / "+r"$\AA$")
axes[0].set_ylabel("In Plane Distance / "+r"$\AA$")
axes[0].set_xticks([0,25,50,75])
#cbar0.set_label('Relative Free Energy / kcal/mol', rotation=270)

pltt=axes[1].contourf(X,Y,free_energy_msm,levels=nlevl,cmap='nipy_spectral')
axes[1].set_xlabel("Out of Plane Distance / "+r"$\AA$")
axes[1].set_xticks([0,25,50,75])
cbar1=plt.colorbar(pltt,ax=axes[1])
cbar1.ax.get_yaxis().labelpad = 15
#cbar1.set_label('Relative Free Energy / kcal/mol', rotation=270)
#plt.scatter(clust_cent_1[:,0],clust_cent_1[:,1],color='k',alpha=0.5)
#plt.xlabel("Out of Plane Distance / Ang")
#plt.ylabel("In Plane Distance / Ang")

pltt=axes[2].contourf(X,Y,free_energy_out,levels=nlevl,cmap='nipy_spectral')
axes[2].set_xlabel("Out of Plane Distance / "+r"$\AA$")
axes[2].set_xticks([0,25,50,75])
cbar2=plt.colorbar(pltt,ax=axes[2])
cbar2.ax.get_yaxis().labelpad = 15
cbar2.set_label('Relative Free Energy / kcal/mol', rotation=270)
plt.tight_layout()

plt.savefig(out_directory+"Figure_3_FES.png",transparent=True,dpi=300)
plt.show()

#exit()

#better to define 3 regions frist
#A: x<=37.0 and any Y
#A_B:x>37.0 and x<=50.0 and any Y
#B: x>=50.0 and x <=60.0 and any Y
#B_C: x>=6 and X <=70 and Y> 5 and Y <22.5
#C: union (x>=60 and X <=70 Y <= 5),(x>=60 and X <=70 Y >22.5), (X>70)
AB_bounds=[45,50]
BC_bounds=[58,63]
region_A=np.argwhere(X<=AB_bounds[0])
region_AB=np.argwhere(np.logical_and(X>AB_bounds[0],X<=AB_bounds[1]))
region_B=np.argwhere(np.logical_and(X>AB_bounds[1],X<=BC_bounds[0]))
region_BC=np.argwhere(np.logical_and(X>BC_bounds[0],X<=BC_bounds[1]))
region_C=np.argwhere(X>BC_bounds[1])

conf_energy=10.0
#find locations where free_energy of each is <=conf_energy kcal/mol (we are using this as regions we are confident, like A in the cartoon)

#shift free_energy_in to miminze difference between FES in this region
def region_error(fes_1,fes_2):
    return np.average(np.abs(fes_1-fes_2))
nshift=1000#~numerical precision on search, more shifts means more accurate
shifts=np.linspace(-10,10,nshift)



locs_msm_in_overlap=region_AB#np.argwhere(np.logical_and(free_energy_msm<=conf_energy,free_energy_in<=conf_energy))
free_energy_msm_in_region=np.array([free_energy_msm[x[0],x[1]] for x in locs_msm_in_overlap if all([free_energy_msm[x[0],x[1]] < conf_energy,free_energy_in[x[0],x[1]] < conf_energy])])
free_energy_in_msm_region=np.array([free_energy_in[x[0],x[1]] for x in locs_msm_in_overlap if all([free_energy_msm[x[0],x[1]] < conf_energy,free_energy_in[x[0],x[1]] < conf_energy])])
#now find locations where msm and out overlap
locs_msm_out_overlap=region_BC#np.argwhere(np.logical_and(free_energy_msm<=conf_energy,free_energy_out<=conf_energy))
free_energy_msm_out_region=np.array([free_energy_msm[x[0],x[1]] for x in locs_msm_out_overlap  if all([free_energy_msm[x[0],x[1]] < conf_energy,free_energy_out[x[0],x[1]] < conf_energy])])
free_energy_out_msm_region=np.array([free_energy_out[x[0],x[1]] for x in locs_msm_out_overlap  if all([free_energy_msm[x[0],x[1]] < conf_energy,free_energy_out[x[0],x[1]] < conf_energy])])

errors_in_msm=[region_error(free_energy_msm_in_region+shift,free_energy_in_msm_region) for shift in shifts]
in_msm_shift=shifts[np.argmin(errors_in_msm)]

#evaluate error between the shifted msm curve (include in_msm_shift) and the shifted out curve
errors_out_msm=[region_error(free_energy_msm_out_region+in_msm_shift,free_energy_out_msm_region+shift) for shift in shifts]
out_in_shift=shifts[np.argmin(errors_out_msm)]

free_energy_msm+=in_msm_shift
free_energy_out+=out_in_shift
#use lowest energy estimate for each at each grid point
stacked_energy=np.stack([free_energy_in,free_energy_msm,free_energy_out],axis=2)
free_energy_final=np.min(stacked_energy,axis=2)


##plot errors to verify optimizaiton
plt.plot(shifts,errors_in_msm,label="msm  in")
plt.plot(shifts,errors_out_msm,label="out  in")
plt.xlabel("Energy Shift from in model / kcal/mol")
plt.ylabel("Overlap Error / kcal/mol")
plt.legend()
plt.tight_layout()
plt.savefig(out_directory+"Figure_Energy_Shift.png",transparent=True,dpi=300)
plt.show()
#exit()


#now use fes_in to start
#free_energy_final=np.zeros_like(free_energy_msm)
#for ind in region_A:
#    free_energy_final[ind[0],ind[1]]=free_energy_in[ind[0],ind[1]]
#for ind in region_AB:
#    free_energy_final[ind[0],ind[1]]=0.5*(free_energy_msm[ind[0],ind[1]]+in_msm_shift + free_energy_in[ind[0],ind[1]])
#for ind in region_B:
#    free_energy_final[ind[0],ind[1]]=free_energy_msm[ind[0],ind[1]]+in_msm_shift
#for ind in region_BC:
#    free_energy_final[ind[0],ind[1]]=0.5*(free_energy_msm[ind[0],ind[1]]+in_msm_shift + free_energy_out[ind[0],ind[1]]+out_in_shift)
#for ind in region_C:
#    free_energy_final[ind[0],ind[1]]=free_energy_out[ind[0],ind[1]]+out_in_shift

free_energy_final=free_energy_final-np.min(free_energy_final)

fig=plt.figure(figsize=(6.5,3.5))
free_energy_final[free_energy_final>=10]=10
plt.contourf(X,Y,free_energy_final,levels=nlevl,cmap='nipy_spectral')
cbar=plt.colorbar()
cbar.ax.get_yaxis().labelpad = 15
cbar.set_label('Relative Free Energy / kcal/mol', rotation=270)
plt.tick_params(axis='both',direction='in')
plt.xlabel("Out of Plane Distance / "+r"$\AA$")
plt.ylabel("In Plane Distance / "+r"$\AA$")
plt.tight_layout()
plt.savefig(out_directory+"Figure_FES.png",transparent=True,dpi=300)
plt.show()

#

msm_in_states=[x[0] for x in np.argwhere(clust_cent_1[msm_1.active_set,0]<=50.)]
msm_out_states=[x[0] for x in np.argwhere(clust_cent_1[msm_1.active_set,0]>=65.)]
print(msm_in_states)
print(msm_out_states)


if recalc_mfpt:
    in_out_mfpt=msm_1.sample_mean('mfpt', msm_in_states, msm_out_states)/convert2
    out_in_mfpt=msm_1.sample_mean('mfpt', msm_out_states, msm_in_states)/convert2
    in_out_mfpt_err=msm_1.sample_std('mfpt', msm_in_states, msm_out_states)/convert2
    out_in_mfpt_err=msm_1.sample_std('mfpt', msm_out_states, msm_in_states)/convert2
    np.save(out_directory+'mfpts_lag'+str(lag)+'_samples'+str(nsampl)+'.npy',[in_out_mfpt,out_in_mfpt,in_out_mfpt_err,out_in_mfpt_err])
else:
    stuff=np.load(out_directory+'mfpts_lag'+str(lag)+'_samples'+str(nsampl)+'.npy')
    in_out_mfpt,out_in_mfpt,in_out_mfpt_err,out_in_mfpt_err=stuff

print("in -> out:",in_out_mfpt,"+/-",in_out_mfpt_err,"ns")
print("out -> in:",out_in_mfpt,"+/-",out_in_mfpt_err,"ns")
