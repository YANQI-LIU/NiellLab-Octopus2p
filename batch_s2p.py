''' Batch run suite2p on acquisitions. Folder selected must be structured as individual folders containing one sbx file and their .mat file each. 
Run in suite2p environment + tkinter module
'''

import os
from pathlib import Path

import suite2p
from suite2p.run_s2p import run_s2p

import tkinter.filedialog
import tkinter as tk

ops = suite2p.default_ops() # populates ops with the default options

ops['batch_size'] = 200 # we will decrease the batch_size in case low RAM on computer

ops['fs'] = 10 # sampling rate of recording, determines binning for cell detection

#ops['tau'] = 0.5 # timescale of calcium indicator to use for deconvolution

ops['input_format'] = "sbx"
# specify that our file input format is sbx

ops["reg_tif"]= True
# save the registreation result as tiff

ops['sbx_ndeadcols']= 70
ops['sbx_ndeadrows']= 4
# an estimate by me, to be refined in the future

#ops['nchannels']=2
# comment this out if only have 1 channels
#ops["reg_tif_chan2"]= True
# save the registered channel 2 tiff, comment out if only 1 channel

# select folder 
data_dir= tk.filedialog.askdirectory(title= 'Select the folder containing individual folders of acquisition' )

db=[]
for i in os.listdir(data_dir):
    db.append({f'data_path': [os.path.join(data_dir,i)]})

for i in db:
    print(i)

user_input= input('Enter 1 to continue:')
if user_input==1:
    for dbi in db:
        output_ops = suite2p.run_s2p(ops=ops, db=dbi)
    print(f'Finished running {len(db)} files!')
    
else:
    exit()
