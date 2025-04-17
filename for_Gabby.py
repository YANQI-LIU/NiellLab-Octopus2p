import os
from pathlib import Path

#import custom library by Anne
from analysis_module import util, constants
import analysis_module as am

# import analysis & plotting related modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from matplotlib.backends.backend_pdf import PdfPages # import pdfpage for saving figures as multipage pdf
import scipy
import pandas as pd
from tqdm import tqdm
from sklearn.preprocessing import minmax_scale
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster 

# import user interface module
import tkinter.filedialog
import tkinter as tk

# define parameters
stim_frames= [199,399,599] # stimulus start on these frames

# Load in s2p result files
#folder= r'D:\Anne\Data_from_Angelique\Data\120123\9\suite2p\plane0\ops0\\'
folder= tk.filedialog.askdirectory(title= 'Select the folder containing s2p results, typically plane\ops0\\' )

F = np.load(os.path.join(folder, 'F.npy'), allow_pickle=True)
Fneu = np.load(os.path.join(folder, 'Fneu.npy'), allow_pickle=True)
#spks = np.load(os.path.join(folder, 'spks_original.npy'), allow_pickle=True)
stat = np.load(os.path.join(folder, 'stat.npy'), allow_pickle=True)
ops =  np.load(os.path.join(folder, 'ops.npy'), allow_pickle=True)
ops = ops.item()
iscell = np.load(os.path.join(folder, 'iscell.npy'), allow_pickle=True)

print('Finished loading acquisition and s2p results')

# Grab master folder, ie. 3 folders up from the ops folder as the output folder
p= Path(folder).parents[2].__str__()
# Instantiating pdf document
PDF = PdfPages(f'{p}\\results.pdf')

# preprocess s2p roi traces (iscell probability filter, normalization of dye loading using median)
fig, fig2= am.plot_s2proi(stat,iscell,ops)


PDF.savefig(fig) # save in pdf
PDF.savefig(fig2)

iscell_thresh= float( input('Select threshold based on s2p cell probability (0 < x <= 1):'))

while not 0 < iscell_thresh <=1:
    print('Must be an integer or float between 0 and 1, please enter again ')
    iscell_thresh= float( input('Define is cell threshold (ie. 0.5):'))

print(f'Selected iscell threshold is {iscell_thresh}')

prob_cell=iscell.T[1]
selected_index= np.argwhere(prob_cell>iscell_thresh)
print(f'There are {len(selected_index)} rois above the probability threshold')

celldff, celldff1= am.norm_fluorescence(F, selected_index)
nroi= celldff.shape[0] # total numbers of rois

fig = plt.figure(figsize=(10, 5))
plt.imshow(celldff, vmin=0, vmax=0.3, aspect='auto',interpolation='none')
plt.title (f'df/f trace (normalized by mediam) for all {celldff.shape[0]} rois identified as cell across all frames')
plt.show(block=False)
plt.pause(0.01)

# save in pdf
PDF.savefig(fig)

# clustering based on the normalized dff
Z, fig= am.response_cluster(celldff)
# hierarchical clustering

# save in pdf
PDF.savefig(fig, bbox_inches='tight')

# ask user input to cut dendrogram, ask again if not happy with cluster number
while True: 
    big_cluster,t = am.define_cluster(Z)
    cluster_confirm= input('Press enter to continue, press anything else to to redefine cutoff:')
    if cluster_confirm =="":
        break
    else:
        print("retry")

# plot the dendrogram with cut line
fig, ax=plt.subplots(figsize=(15, 5))
dn = dendrogram(Z)
ax.axhline(t,ls= '--', c='k')
# plot a horizontal cut off line to the dendrogram

# save in pdf
PDF.savefig(fig)
plt.show(block=False)
plt.pause(0.01)

# plot roi on mean image color coded by cluster
cellstat=stat[selected_index.flatten()]
# stat from structures pre-selected (ie. by a defined threshold of iscell probability)

cellstat=cellstat[~np.isnan(celldff1.min(axis=1))]
# also remove those indices that were dividing by zero during our celldff process

fig= am.plot_clust_roi(big_cluster, cellstat,ops)
# plot out the roi location

# save in pdf
PDF.savefig(fig)


##### detailed description of clusters
for i in range(1,big_cluster.max()+1):

    # prepare data for plotting
    cluster_cells= cellstat[big_cluster==i]
    cluster_response= celldff[big_cluster==i]

    # create subplots
    fig = plt.figure(figsize=(24,8))
    ax1= plt.subplot(1,4,1)
    ax2= plt.subplot(1,4,(2,3))
    ax3= plt.subplot(1,4,4)

    center_row= []
    center_col = []
    
    for j in cluster_cells:
        center_row.append( j['med'][0])
        center_col.append( j['med'][1])
    
    ax1.imshow(cluster_response, aspect='auto', vmin=-1, vmax=1, interpolation='None')
    ax1.set_xlabel('frame #')
    ax1.set_ylabel('nroi')
    ax1.set_title("roi dff within this cluster")

    
    ax3.imshow(ops['meanImg'], cmap='gray')
    ax3.scatter(center_col, center_row, s= 3, alpha= 0.8, color='r')
    ax3.set_title("roi location within this cluster")
    ax3.axis('off')

    # prepare for the next page- staggered roi dff

    # initialize empty array to map staggered responses
    stagger_df= np.zeros_like(cluster_response)
    for k in range(0,cluster_response.shape[0]):
        this_roi= cluster_response[k,:]+ k*0.1 
        #add an arbitrary value so the traces are not on top of each other
        stagger_df[k,:]= this_roi

    
    ax2.plot(stagger_df.T) # need to plot the transposed array because plt.plot consideres columns to be lines
    ax2.vlines(x=stim_frames,ymin=-0.1 ,ymax=(k+1)*0.1, color='k',linestyle='--')
    ax2.set(yticklabels=[])
    ax2.set_xticks(stim_frames, ['Stim1','Stim2', 'Stim3'])
    ax2.set_title("roi dff within this cluster")

    fig.suptitle(f'Cluster # {i}')

    PDF.savefig(fig) # save in pdf INSIDE for loop so each for loop is one page!!!
    
    plt.show(block=False)
    plt.pause(0.01)

# plot distribution % and numbers of each ROI in clusters
unique, counts = np.unique(big_cluster, return_counts=True)
fig= plt.figure(figsize=(8, 8))
plt.pie(counts, labels=unique, autopct=lambda p:f'{p:.1f}% ({p*sum(counts)/100 :.0f}) ')
plt.title('ROI distribution in each cluster')

PDF.savefig(fig)
plt.show(block=False)
plt.pause(0.01)

# close the PDF
PDF.close()