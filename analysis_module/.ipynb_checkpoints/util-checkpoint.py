''' useful functions and constants for analysis of 2p ROIs!

'''

import numpy as np
import pandas as pd
from matplotlib import colormaps as cm


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
        
# create object structure

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
    