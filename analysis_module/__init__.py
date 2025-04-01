''' useful functions for analysis of 2p ROIs!

'''

import os

# import sbx modules
from sbxreader import sbx_get_metadata
from sbxreader import sbx_memmap
from sbxreader import sbx_get_info

import tkinter.filedialog
import tkinter as tk
# import user interface module

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster 

from tqdm import tqdm
import seaborn as sns


def get_sbx(file_name):

    ''' read sbx metadata
    '''
    metadata = sbx_get_metadata(file_name.name)
    #print(metadata) # dictionary with the recording metadata
    
    info = sbx_get_info(file_name.name)
    # load mat file using scipy.io as a matlab_structure
    
    return metadata, info

def parse_events(stimulus): 
    '''  get time points of each stimulus and re-link back to imaging frames
    almost exact copy pasta from get2pdata_sbx.m, including the comments '''
    
    TTL0 = 1
    TTL1 = 2
    TTLboth =3
    # new TTL format, contains both rising and falling
    
    if max(stimulus.event_id)> 3:
        stimulus.event_id[stimulus.event_id==1] = TTL0 # TTL0 rising
        stimulus.event_id[(stimulus.event_id==4)|(stimulus.event_id==8)] = TTL1 # TTL1 rising or falling
        stimulus.event_id[(stimulus.event_id==5)|(stimulus.event_id==9)] = TTLboth # TTL0 rising  + TTL1 rising
    
    stim= np.where((stimulus.event_id==TTL0)|(stimulus.event_id==TTLboth))
    fr= (stimulus.frame[stim] + stimulus.line[stim]/796)
    # get frame each trigger occured on, add fraction based on line
    # ??? not sure if used later
    
    # get phase trigger signals
    phasesync= np.where((stimulus.event_id==TTL1)|(stimulus.event_id==TTLboth))
    if stimulus.frame[phasesync[0][0]]==0:
        phasesync= phasesync[1:]
        # something funny happens at start giving transition on the first frame
    phasesync=phasesync[0][0::2]  #scanbox records rising and falling edge;
    fr1= (stimulus.frame[phasesync] + stimulus.line[phasesync]/796) #get frame each trigger occured on,  add on fractio nbased on line
    phasetimes= fr1*stimulus.dt
    # this is the actual time where each stimulus starts in seconds
    
    # get phase trigger signals
    frsync= np.where((stimulus.event_id==TTL0)|(stimulus.event_id==TTLboth))
    fr2= (stimulus.frame[frsync] + stimulus.line[frsync]/796)
    videoframetimes= fr2*stimulus.dt
    # videoframetimes used for spars noise stimulus
    return fr, phasesync, phasetimes, frsync, videoframetimes

def parse_stimulus(stimulus, phasetimes): 
    ''' Parse out events based on the recording.mat file
    '''
    # For non-sparse noise stimulus
    cycLength=np.median(np.diff(phasetimes))/stimulus.dt # number of frames in window around each cycle, min of 4 sec, or actual cycle lengt+2
    # so this is the number of frames around each stimulus 
    
    cycWindow= int(np.round(np.max([2/stimulus.dt, cycLength])))
    # get the higher value between cycle length (calculated based on 2 seconds), or actual cycle length
    
    startFrame=int(np.round( (phasetimes[0]-1)/stimulus.dt)) # movie starts 1 second before first stim
    stimTimes= np.array(np.arange(1,(stimulus.nframes-startFrame-cycWindow-30)*stimulus.dt,cycLength*stimulus.dt)) 
    #make sure you have one cycle, plus 3sec padding to be safe, at end
    
    stimTimesOld = phasetimes-phasetimes[0]+1
    stimTimes=stimTimes[stimTimes<np.max(stimTimesOld)]
    # Not sure what these two steps does?
    
    ncycles= len(stimTimes)
    # total numbers of individual stimulus presented
        
    return cycWindow,cycLength, startFrame, ncycles,stimTimes

def parse_sn(videoframetimes, stimulus):

    ''' parse out stimulus times for sparse noise stimulus

    '''

    cycLength=np.median(np.diff(videoframetimes))/stimulus.dt # number of frames in window around each cycle, min of 4 sec, or actual cycle lengt+2
    # so this is the number of frames around each stimulus 
    cycwindow= int(np.round(np.max([2/stimulus.dt, cycLength])))
    # get the higher value between cycle length (calculated based on 2 seconds), or actual cycle length

    startFrame=int(np.round( (videoframetimes[0]-1)/stimulus.dt)) # movie starts 1 second before first stim
    stimTimes= videoframetimes-videoframetimes[0]+1

    stimTimes= stimTimes[stimTimes<(stimulus.nframes-startFrame)*stimulus.dt-5] 
    stimFrames= np.round(stimTimes/stimulus.dt)

    framerange=np.arange(0,1+cycwindow*2)

    return cycwindow, startFrame, stimFrames, framerange

def get_stimorder(stimulus,acquisition, cycLength, ncycles):
    ''' get the order of stimulus from acquicision info and cycle length
    '''
    
    # to fill in the zeros between stimulus as the previous stimulus. especially important step to process 6x4 on/off stimulus
    stimRec_cond=[]
    # initiate empty list to store
    
    for count, item in enumerate(acquisition['stimRec']['cond'][0][0]):
        if item ==0:
            stimRec_cond.append(acquisition['stimRec']['cond'][0][0][count-1][0])
        else:
            stimRec_cond.append(item[0])
    # so annoying when .mat structure files gets loaded as lits within list within list!!
    
    # process acquisition session data
    
    # the way that .mat file are read is very annoying- its arrays within arrays within arrays
    stimT=acquisition['stimRec']['ts'][0][0]- acquisition['stimRec']['ts'][0][0][0]
    # so variables in acquisition['stimRec'] are in units of frames on the stimulus display (which is x seconds* 60Hz)
    
    # generate a list of stimulus order, found by comparing where the stimulus falls at the begining of each cycle in seconds
    stimorder=[]
    for i in range(0,ncycles):
        stim_index= np.min(np.where(stimT>((i)*cycLength*stimulus.dt+0.1))[0])
        stimorder.append(stimRec_cond[stim_index])
    
    nstim=np.max(stimorder)
    print(f'There are {nstim} different stimulus')
    
    return stimorder

def get_stim_frame(stimulus, stimorder, cycLength, startFrame, stimTimes, cycWindow):
    ''' assign a stimulus event to each frame of the recording
    '''

    # get each stimulus id and how many times its been repeated 
    nstimrep=[]
    # initiate a list to store reptition of stimulus
    for i in np.unique(stimorder):
        this_stim_rep= np.sum(stimorder==i)
        nstimrep.append(this_stim_rep)
        # fill in the list, each stimulus (index) is repeated x times
    
    total_frames= cycLength*np.max(stimorder)
    reps= np.floor((stimulus.nframes-startFrame)/total_frames)
    print(f'Full stimulation set repeated for {reps} times')

    # assign a stimulus event to each frame
    stimframe= np.zeros((stimulus.nframes))
    # initiate a array of zeros to a stimulus event for each frame
    
    stim_id= []
    stim_rep=[]
    stim_start=[]
    stim_end=[]
    # initiate array to store stimulus identifier and repetition window
    
    for count,element in enumerate(np.unique(stimorder)):
        replist=np.where(stimorder==element)
        # find the corresponding repetition times for this stimulus
        
        for j in np.arange(0, nstimrep[count]):
            # for each repetition, gather the start time, fill in stimulus identifier for the duration of cycle
            start=stimTimes[replist[0][j]]/stimulus.dt
            window_start= np.round(start)+startFrame # dont for get that there are startFrames number of frames before the stimulus actually starts!
            window_end= window_start+ cycWindow-1
            #this_window= np.round(start + np.arange(1,cycwindow))
            stimframe[int(window_start):int(window_end)]=element
            
            stim_id.append(element)
            stim_rep.append(j)
            stim_start.append(int(window_start))
            stim_end.append(int(window_end))

    #store results in pandas dataframe
    d={'stim_id':stim_id, 'stim_rep':stim_rep, 'stim_start':stim_start, 'stim_end':stim_end}
    stim_pd= pd.DataFrame(d)
    
    return stim_pd, nstimrep

def plot_s2proi(cellstat,iscell, ops):

    ''' plots all roi on the mean image colorcoded based on their is_cell probability, and plot roi with prob intervals on mean image
    for setting a threshold to include in the subsequent analysis
    Return the two figures produced for saving in pdf
    '''

        # plot all roi 
    center_row= []
    center_col = []
    prob_cell=iscell.T[1]
    
    for j in cellstat:
        center_row.append( j['med'][0])
        center_col.append( j['med'][1])
    
    fig = plt.figure(figsize=(6, 6))
    plt.imshow(ops['meanImg'], cmap='gray')
    plt.scatter(center_col,center_row, s= 10, alpha= 0.8, c=prob_cell)
    plt.title('All Rois on mean image color coded by is_cell probability')
    plt.colorbar()

    # plot roi based on their probability intervals
    fig2, axs = plt.subplots(3,3,figsize=(12,12))
    axs=np.ravel(axs)
    
    for count,i in enumerate(np.arange(0.1,1,0.1)):
        this_index= np.argwhere((prob_cell<(i+0.1))&(prob_cell>i))
        this_cellstat=cellstat[this_index]
    
        center_row= []
        center_col = []
        
        for j in this_cellstat:
            center_row.append( j[0]['med'][0])
            center_col.append( j[0]['med'][1])
        
        axs[count].imshow(ops['meanImg'], cmap='gray')
        axs[count].scatter(center_col,center_row, s= 10, alpha= 0.8)
        axs[count].set_title(f' {i:.2f} >=cell prob <{(i+0.1):.2f}')
        axs[count].axis('off')
    
    axs = np.reshape(axs, (3,3))
    plt.subplots_adjust(hspace= 0.1,wspace= 0.1)
    return fig, fig2

def norm_fluorescence(F, selected_index):
    ''' Median filter each roi fluorescent trace followed by normalizing (by median, to control for dye loading variations)
    '''
    # median filter fluorescent trace for all rois
    F_filtered= scipy.signal.medfilt(F, [1,3])
    # median filter of size 3 frame along the 2nd dimension, which is the time series for each roi
    
    cellF=F_filtered[selected_index.flatten()]
    # take cells defined by iscell threshold from previous section
    
    print(f'there are {cellF.shape[0]} cells out of total identified structure {F.shape[0]}')
    
    np.seterr(divide='ignore', invalid='ignore')
    # somehow there are cells with fluorescent values = nan...??? take cares of division by zero cases..
    
    baseline= np.median(cellF, axis=1) # use median value of each roi trace as baseline
    celldff1= (cellF - baseline[:, np.newaxis])/baseline[:, np.newaxis]
    # calculate DFF change in fluorescence using median value (to control for amount of dye loading)
    
    celldff=celldff1[~np.isnan(celldff1.min(axis=1))]
    # remove traces that was divided by zero

    print(f'excluding weird fluorescent traces that have nans as median, the remaining good traces are {celldff.shape[0]}')
    
    #here we also store celldff1 for later use to exclude rois with median values=0 
    return celldff, celldff1

def reorganize_dff(celldff, cycwindow, nstim, nstimrep, stim_pd, offset, baseRange):
    ''' organize traces for each roi in to each cycle for each stimulusid and for each repetition,
    then subtract baseline mean from each trace
    '''

    reorganized_trace= np.zeros((celldff.shape[0],cycwindow-1, nstim, np.max(nstimrep)))+np.nan
    # create storage to store mean images at a given stimulus at a given repetition
    print(f'cell dff trace array reorganized to {reorganized_trace.shape}, with ncells, cycwindow-1, stim_id, max repetition')

    # extract evoked responses based on pre-defined periods
    for count, element in enumerate(tqdm(np.unique(stim_pd['stim_id']))):
        for j in range(0, nstimrep[count]):
            # or np.arange(0, int(reps)) when not all stim id is repeated the same time (ie some have 1 less)- NEED TO THINK ABOUT THIS
            
            this_pd= stim_pd[(stim_pd['stim_id']==element) & (stim_pd['stim_rep']==j)]
            this_trace= celldff[:, this_pd['stim_start'].values[0]+offset : this_pd['stim_end'].values[0]+ offset]
            # grab relevant frames for a given stimulus at a given repetition

            this_baseline= np.nanmean(this_trace[:, baseRange], axis=1)

            trace_evoked= this_trace - this_baseline[:,np.newaxis]
            # calculate evoked image by using evoked frame ranges and baseline frame ranges, then mean this these frames
            
            reorganized_trace[:,:,count,j]= trace_evoked
            # store into storage
            
    return reorganized_trace

def plot_avgtrace(reorganized_trace, color_dict, style_dict,alpha_dict):

    avg_stim_trace= np.nanmean(reorganized_trace, axis= 3)
    # trace for each stimulus- averaged across all repititions

    fig=plt.figure()

    for i in np.arange(0,avg_stim_trace.shape[2]):
        plt.plot(np.mean(avg_stim_trace[:,:,i], axis=0), label=f'stimulus # {i}',color= color_dict[i], linestyle= style_dict[i], alpha= alpha_dict[i])
    plt.title('mean response across all rois for each stimulus')
    plt.legend(loc='best',bbox_to_anchor=(1.1, 1.05))
    plt.show()
    return fig

def find_responsive(trace, base_range,evo_range, nstim, tail='greater'):

    ''' find rois that are significantly responsive to a certian stimuli (takes in acount for all repetitions for a given stimulus) 
    inputs dfof trace that is rearranged as (ncells, cycwindow-1, stim_id, repitition), base frame range and evoked frame range, and tail ('greater' or 'less')
    '''

    response_matrix= np.zeros((trace.shape[0],trace.shape[2]))
    # create empty array with shape (nroi, nstim)
    
    for i in tqdm(range(0,trace.shape[2])):
    # loop over different stimulus
    
        for j in range(0,trace.shape[0]):
        # loop over different rois
            bl_measures= trace[j,base_range,i,:]
            evo_measures= trace[j,evo_range,i,:]
            ranksums= scipy.stats.ranksums (evo_measures.flatten(),bl_measures.flatten(), alternative=tail)
    
        # record wether this roi during evoked response period is significantly greater than during baseline period determined by the ranksums test
            if ranksums.pvalue <(0.05/nstim): # Bonferroni correction to control for multiple comparisons
            #if compare large enough numbers of features then there is morelikely that something will appear difference by chance
                response_matrix[j,i]=1 # assign those fall below p value as 1
        else:
            pass
    
    return response_matrix

def response_cluster(cluster_input_responsive):
    ''' hierarchical clustering based on response from statistial test
    '''

    Z = linkage(cluster_input_responsive, 'ward')
    plt.figure(1)
    plt.subplots(figsize=(15, 5))
    dn = dendrogram(Z)
    plt.show()

    
    sns_figure=sns.clustermap(
    cluster_input_responsive,
    figsize=(5,5),
    method='ward',
    row_cluster=True,
    col_cluster=False,
    cmap='viridis')
    plt.title('whether each unit is significant to a stimulus or not')
    plt.show()

    return Z, sns_figure.figure

def response_cluster_sn(ndarray):
    ''' hierarchical clustering based on response from statistial test
    '''

    Z = linkage(ndarray, 'ward')
    plt.figure(1)
    plt.subplots(figsize=(15, 5))
    dn = dendrogram(Z)
    plt.show()

    
    sns_figure=sns.clustermap(
    ndarray,
    figsize=(5,5),
    method='ward',
    row_cluster=True,
    col_cluster=False,
    cmap='viridis',
    vmin = 0,
    vmax= 0.3)
    plt.title('mean dff during taurange in each frame for each roi')
    plt.show()

    return Z, sns_figure.figure

def define_cluster(Z):
    t= float( input('Input cut off level based on dendrogram (int or decimal >0): '))

    big_cluster= fcluster(Z, t=t, criterion='distance')
    print(f'number of arbitrary cluster= {big_cluster.max()}')

    return big_cluster,t

def plot_clust_roi(big_cluster, cellstat, ops):

    ''' Plot location of roi based on cluster result
    '''
    print("Plotting positions of rois in each cluster")
    fig=plt.figure()
    
    for i in tqdm(range(1,big_cluster.max()+1)):
        cluster_cells= cellstat[big_cluster==i]
        
        center_row= []
        center_col = []
        c=['pink', 'olive', 'cyan','blue', 'orange', 'green', 'red', 'purple', 'brown']
        
        for j in cluster_cells:
            center_row.append( j['med'][0])
            center_col.append( j['med'][1])
            
            plt.imshow(ops['meanImg'], cmap='gray')
            plt.scatter(center_col, center_row, s= 5, alpha= 0.8, c= c[i-1])
    plt.title('Clustering based on stimulus significance result')
    plt.show()
    return fig
        
