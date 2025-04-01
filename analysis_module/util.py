''' useful functions for analysis of 2p ROIs!
Includes functions specific to each stimulus

'''

import numpy as np
import pandas as pd
from matplotlib import colormaps as cm

import scipy

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

# 8 way stimulus parameters and graphic styles
def gratings_param(acquisition, stimorder):
    ''' Load stimulus parameters and set graphic styles for plotting

    '''
    print('Current stimulus type: 8 ways grating stimulus, total numbers of stimulus: 17')

    orient=acquisition['orient']
    orient=np.append(orient, np.nan)
    # the last stimulus for this type of stimulus session is a full field flash, so append nan for orientation
    
    freq=acquisition['freq']
    freq=np.append(freq, np.nan)
    # the last stimulus for this type of stimulus session is a full field flash, so append nan for spatial frequency

    d= {'stim_id': np.unique(stimorder), 'orientation':orient, 'spatial_frequency':freq}
    stim_detail= pd.DataFrame(d)

    # create velocity vectors - reduces the orientation dimension but adds a third dimension
    direction_dict= {'orientation': [0, 45,90,135,180,225,270,315],
                    'direction_x': [1,1,0,-1,-1,-1,0,1,],
                    'direction_y': [0,1,1,1,0,-1,-1,-1]
                    }
    direction_pd= pd.DataFrame(data=direction_dict) 
    stim_detail= pd.merge(stim_detail, direction_pd, on='orientation', how='left') # merge data frame 

    cmap=cm['Accent'].colors

    color_dict= np.repeat(cmap,2,axis=0)
    color_dict=np.append(color_dict,[[0.1,0.1,0.1]], axis=0)
    # assign colors to stimulus- stimulus with same orientation gets same color, black as full field flash

    style_dict= ['dashed','dotted','dashed','dotted','dashed','dotted','dashed','dotted','dashed','dotted','dashed','dotted','dashed',
                 'dotted','dashed','dotted','solid']
    # assign line style as indications of spacial frequency- dashed correspond to sf= 0.01, dotted sf=0.16, solid line as full field flash

    return stim_detail, color_dict, style_dict


def load_snmovie():
    ''' Load and preprosess sparse noise movie
    returns modified m (with full field flash removed), m_full (array of original movie) and size_m (array indicating size)
    '''

    movie= scipy.io.loadmat(r'D:\Anne\Matlab_test\data\octo_sparse_flash_10min.mat') # load movie data
    m=(movie['moviedata'].astype('float')-127)/128  # convert to float (correspond to double() in matlab),  normalize pixel values 
    m_full=m.copy() # full version of the movie

    m[movie['sz_mov']==255]=0 # something about removing full field..? copied from matlab code. 
    #Probably have something to do for calculating STA so full field flashes doesn't compromise the appearance of STA (ie. too bright)
    sz_m= movie['sz_mov']
    # movie['sz_mov'] is an array with the same size as movie with value at each pixel indicate size of connected component

    return m, m_full, sz_m

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
    temp=cluster_df.groupby(by='orientation').sum().reset_index() # sum up weights that belong to the same orientation regardless of spatial frequency
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
    temp=cluster_df.groupby(by='orientation').sum().reset_index() # sum up weights that belong to the same orientation regardless of spatial frequency
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

def percentile_sta(resp, sta, percentile, nframes, high=True):
    ''' calculate mean sta from frames with responses that are above or below x percentile
    stimulus triggered average, related to sparse noise like stimulus
    '''
    if high== True:
        this_index= np.argwhere(resp[:nframes]>=np.percentile(resp, percentile))
    else:
        this_index= np.argwhere(resp[:nframes]<=np.percentile(resp, percentile))
    # find time points where response is greater than 75% percentile or less than 25% percentile
    # crop the length to nframe (in this case video timepoints are shorter than recording time points)
    
    this_sta= sta[:,:,this_index.flatten()]
    this_mean_sta= np.mean(this_sta,2)
    return this_mean_sta 

def find_frame_sz(sz_m, m_full, max_ind, nframes):
    ''' find frames at a given pixel location of varying size and intensity
    related to sparse noise like stimulus
    '''
    
    sizes= np.unique(sz_m) # find all unique sizes
    sizes=sizes[1:] # remove the first value (0)
    
    # crop out non-recorded frames from both intensity move and size movie
    m_full=m_full[:,:,:nframes]
    sz_m= sz_m[:,:,:nframes]
    
    size_frame_on=[] # initialize list to store frames
    for i in sizes:
        these_frame= np.where((sz_m[max_ind[0], max_ind[1],:]== i) & (m_full[max_ind[0], max_ind[1],:]== 1) ) # find a given size for ON sitmulus (ie. pixel value in m is 1)
        size_frame_on.append(these_frame)
        
    size_frame_off=[] # initialize list to store frames
    for i in sizes:
        these_frame= np.where((sz_m[max_ind[0], max_ind[1],:]== i) & (m_full[max_ind[0], max_ind[1],:]< 0) ) # find a given size for ON sitmulus (ie. pixel value in m is 1)
        size_frame_off.append(these_frame)
        
    return size_frame_on, size_frame_off

def plot_sz_effects(this_cycle, size_frames, axs, intensity ):
    ''' plot size effects of mean responses across roi
    each transparent dashed line represents a frame with a given size and intensity
    solid line represents mean response across frames

    related to sparse noise like stimulus
    '''

    cmap= ['r','g','b','k']
    labels= ['2', '4', '8', 'full']
    
    for count, item in enumerate(size_frames):
        # so far not including the full field flash (to be fixed later) 
        
        this_size_cycle= this_cycle[item,:,].squeeze()
        mean_cycresp= np.mean(this_size_cycle,axis=2)
        # mean across rois, leaving a nframe x frame range array
        axs.plot(mean_cycresp.T, color= cmap[count], linestyle= '--', linewidth=0.5, alpha= 0.5) # need to plot the transposed array
        axs.plot(np.mean(mean_cycresp, axis=0),color= cmap[count], label=labels[count])
        axs.legend()
        axs.set_title(f'Size effect -  {intensity}')
    
    return 