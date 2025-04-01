''' useful constants for analysis of 2p ROIs!
directly taken from matlab version 

'''
import numpy as np

# Define parameters
offset= -10 
# hard coded in matlab version, give 10 frames buffer in front of each cycle

baseRange=range(0,5)
# take the first five frame as the baseline

evRange= range(9,19)
# take the 9th through 19th frame as evoked activity

# define parameters for sparse noise recordings
offset_sn= 0
baseRange_sn=range(0,2)
evRange_sn= range(9,19)

#other parameters for for octo_sparse_flash_10min.mat
tau= 8
tau_range=np.arange(-4,5) # goes from -4 up to 4 (not including 5)
crange=[-0.1,0.1] # not even sure if its useful later but keep it here for now
