''' useful functions for analysis of 2p ROIs!

'''

import numpy as np
import pandas as pd
from matplotlib import colormaps as cm

import seaborn as sns

# create object structure to store sbx file parameters
class Getdata:
    def __init__(self,metadata, info):
        
        self.resfreq= 7920 # defined by cmn 072519 because it is not 8000 as set in the metadata
        self.framerate= 7920/ info.config.lines
        self.dt= np.round(1/self.framerate, decimals=4)
        # info about frame rates and time in seconds between each frame
        
        self.size= info.sz
        self.nframes= metadata['num_frames']
        # info about x, y and z of image stack
        
        self.event_id= metadata['event_id']
        self.frame= info.frame
        self.line= info.line
        # info about event triggers
        

class configuration: # not useful for now
    def __init__(self, spatialBin, temporalBin, tau, tau_range):
        self.spatialBin=spatialBin
        self.temporalBin= temporalBin
        self.tau=tau
        self.tau_range=tau_range


# 6x4 stimulus parameters and graphic styles
def retinotopy_param(acquisition, stimorder):
    ''' Load stimulus parameters and set graphic styles for plotting
    '''
    print('Current stimulus type: 4 x 6 ON/OFF stimulus, total numbers of stimulus: 48')
    
    positionX=acquisition['positionX'].flatten()
    positionY=acquisition['positionY'].flatten()
    contrast=acquisition['contrast'].flatten()
        # where -1 is off and 1 is on..?

    #store stim id, positionx positiony and contrast  in dataframe
    d= {'stim_id': np.unique(stimorder), 'positionX':positionX, 'positionY':positionY, 'contrast':contrast}
    stim_detail= pd.DataFrame(d)

    # define color map and linestyles for plotting
    # assign colors to stimulus- stimulus with different x position vary by color color, stimulus  position y varies by alpha value
    
    cmap=cm['Accent'].colors
    color_dict_temp= np.repeat(cmap[0:6],4,axis=0) # take the first 6 colors (pos X change), and each repeat four times (posY change)
    color_dict= np.concatenate((color_dict_temp, color_dict_temp))
    
    style_dict= ['dashed'] *24 + ['solid']* 24
    # assign line style as indications contrast, with dashed = -1 and solid = 1
    
    alpha_dict=(0.1,0.4,0.7,1)*12
    # four transparency level,  smaller (negative) y position correspond to more transparent, largest y position is not transparent
    
    marker_dict= ['_']*24+['P']*24
    # markder style represents off(-) or on (+) stimulus

    return stim_detail, color_dict, style_dict, alpha_dict, marker_dict

# 6x4 stimulus parameters and graphic styles
def gratings_param(acquisition, stimorder):
    ''' Load stimulus parameters and set graphic styles for plotting
    '''
    return 

# define useful functions for data-preparation and plotting

def permute_cluster_response(cluster_response, npermute, rng_seed):
    ''' 
    Permute an array of cluster_response of shape nrois x response for npermute numbers of time.
    each row is permuted npermute numbers of time and the average is noted in a temporary array
    the final average across averages and std of all nroi is stored as chance_mean and chance_std

    '''
    cluster_response_chance= np.zeros((cluster_response.shape[1], npermute)) # create empty array to hold sample shuffling
    
    for j in range(0, npermute):
        cluster_response_chance[:,j]=np.mean(rng_seed.permutation(cluster_response, axis=1), axis=0)
    # permute sample x times along the nstim axis

    chance_mean=np.mean(cluster_response_chance, axis=1) # compute mean response strength for shuffled data
    chance_std= np.std(cluster_response_chance, axis=1) # compute std response strength for shuffled data
    
    return chance_mean, chance_std

def arrange_cluster_df(cluster_response, chance_mean,chance_std, stim_detail):
    ''' 
    calculate cluster mean, standard deviation, compares cluster mean with chance_mean plus or minus chance_std
    stores cluster mean, cluster std, and 'weights' of important stimulus in a data frame
    'weights' of important stimulus- differences for cluster means that are above and chance_mean+std and those that are below chance_mean-std
    '''
    cluster_mean= np.mean(cluster_response, axis=0)
    cluster_std=np.std(cluster_response, axis=0)
    
    dif_care_response= cluster_mean- (chance_mean + chance_std)
    dif_care_response[dif_care_response<0]=0 # assign zero for all those differences that are smaller than the mean+ std
    #dif_care_response= dif_care_response.reshape(2,nstim) # reshape to 2 x nstim array so the second set is for sig decrease in fluorescence
    
    dif_notcare_response= cluster_mean- (chance_mean - chance_std)
    dif_notcare_response[dif_notcare_response>0]=0 # assign zero for all those differences that are greater than the mean- std
    #dif_notcare_response= dif_notcare_response.reshape(2,nstim)
    # currently dif_dec_response is not being plotted
    
    cluster_df=stim_detail.copy()
    cluster_df['response_mean']= cluster_mean
    cluster_df['response_std']= cluster_std

    
    cluster_df['care_sig_inc_weight']= dif_care_response
    cluster_df['care_sig_dec_weight']= dif_notcare_response
    
    return cluster_df

def plot_amp_ev(celldf, evoked_range, xval, color_dict, axs=None):
    ''' calculate mean amplitude during evoked period of each unit for each stimulus
    Input cell fluorescent traces, arranged as nroi x nframes(in a cycle) x stimulus
    '''
    stim_resp_amp=[]
    stim_resp_sem=[]
    
    for i in range(0, celldf.shape[2]):
        cluster_fevmean=np.mean(celldf[:,evoked_range, i], axis=1)
        # calculate mean fluorescent during the evoked frames for each units in this cluster
        
        stim_resp_amp.append(np.mean(cluster_fevmean))
        stim_resp_sem.append(np.std(cluster_fevmean)/np.sqrt(celldf.shape[0]))
    
    axs.bar(xval, stim_resp_amp, color=color_dict, yerr=stim_resp_sem)

    return


def plot_bar_appearance(cluster_df, col='',axs=None, cmap=cm['Accent'].colors):
    ''' for 8 way gratings
    Plot what the bar look like for each stimulus that is important to the roi
    need to specify which column to use (ie. care_sig_inc_weight) and to plot on which axs
    '''

    # generate some sort of representation for these stimulus
    center = [5, 3]
    L = np.array([-5, 5]) # length of line
    temp=this_clusterpd.groupby(by='orientation').sum().reset_index() # sum up weights that belong to the same orientation regardless of spatial frequency
    #cmap=cm['Accent'].colors

    for i in range(0, len(temp)):
        phi = np.deg2rad(temp.iloc[i].orientation+90) # so that degree 0 is vertical 
        x = center[0] + np.cos(phi) * L
        y = center[1] + np.sin(phi) * L
    
        this_weight= temp.iloc[i][col]
        axs.plot(x, y, linewidth=this_weight*50, color=cmap[i],alpha=0.75,
                    label=f'{temp.iloc[i].orientation}')        
    return

def plot_direction(cluster_df,col='',axs=None, cmap=cm['Accent'].colors):
    ''' for 8 way gratings 
    Plot what the motion direction like for each stimulus that is important to the roi
    need to specify which column to use (ie. care_sig_inc_weight) and to plot on which axs
    '''    
    temp=this_clusterpd.groupby(by='orientation').sum().reset_index() # sum up weights that belong to the same orientation regardless of spatial frequency
    #cmap=cm['Accent'].colors
    for i in range(0, len(temp)):
        
        this_x= temp.iloc[i].direction_x 
        this_y= temp.iloc[i].direction_y 
    
        this_weight= temp.iloc[i][col]
    
        axs.arrow(0,0,this_x,this_y, head_width=np.abs(this_weight), linewidth=this_weight*20, color=cmap[i],
                     label=f'motion_vector= {this_x,this_y},  {temp.iloc[i].orientation} ') 
    
    return



def plot_scatter_kde_all(df,color_dict, axs=None):
    ''' for 6x4
    plots scatter plot on locations- positionX, position Y, with size of marker the weight of response
    approximate a kernel density estimation based on the weights
    positive weight-  red, negative weight blue
    '''
    weight_positive= df['response_mean'].copy()
    weight_positive[weight_positive<0]=0
    
    weight_negative= df['response_mean'].copy()
    weight_negative[weight_positive>0]=0
    
    axs.scatter(x=df["positionX"], y=df["positionY"], s=20, c='k', alpha= 0.2, marker='x')
    # plot all possible locations, using x as markers
    
    axs.scatter(x=df["positionX"], y=df["positionY"], s=weight_positive*500, c=color_dict[0:24], edgecolors='k',marker='^')
    #plot only significant locations, using circles and size represents the mean responsive weight
    sns.kdeplot(data=df, x='positionX', y='positionY', weights=weight_positive, cmap="rocket", levels=5, thresh=0.4, ax=axs)
    
    axs.scatter(x=df["positionX"], y=df["positionY"], s=abs(weight_negative)*500, c=color_dict[0:24], edgecolors='k',marker='v')
    #plot only significant locations, using circles and size represents the mean responsive weight
    sns.kdeplot(data=df, x='positionX', y='positionY', weights=np.abs(weight_negative), cmap="mako", levels=5, thresh=0.4, ax=axs)

    return 

def plot_scatter_sig(df,color_dict, axs=None):
    ''' for 6x4
    plots scatter plot on locations- for those that are significant only!
    - positionX, position Y, with size of marker the weight of response
    positive weight-  red, negative weight blue
    '''
    
    axs.scatter(x=df["positionX"], y=df["positionY"], s=20, c='k', alpha= 0.2, marker='x')
    # plot all possible locations, using x as markers
    
    axs.scatter(x=df["positionX"], y=df["positionY"], s=df['care_sig_inc_weight']*500, c=color_dict[0:24], edgecolors='k', marker='^', label='above')
    #plot only significant locations, using circles and size represents the mean responsive weight
    
    axs.scatter(x=df["positionX"], y=df["positionY"], s=np.abs(df['care_sig_dec_weight'])*500, c=color_dict[0:24], edgecolors='k', marker='v', label='below')
    #plot only significant locations, using circles and size represents the mean responsive weight

    return 
