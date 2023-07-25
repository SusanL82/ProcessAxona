# -*- coding: utf-8 -*-
"""
Reads axona .bin files from multiple sessions,
then merges the files, downsamples to 6kHz, saves a .mat file for each tetrode.
"""

import tqdm as tqdm #I'm not sure why this is even needed, but it stops an error?
import spikeinterface as si
import spikeinterface.extractors as se
from scipy.io import savemat
import numpy as np
from scipy import signal

##
#allnames = ["pre"]
#allsess = [[1,2,3,4,5]]

allnames = ["pre","post1","post2a", "post2b","post2c"]
allsess = [[1,2],[4,5],[7,8,9,10],[11,12,13,14],[15,16,17,18]]

#allnames = ["post2a", "post2b","post2c"]
#allsess = [[13,14,15,16],[17,18,19,20],[21,22,23,24]]
Tetnums = [2,4,5,6,7,8,9,11]

InFolder = "E:/SEPSIS/RAW/BIN/HD214-2004/"
OutFolder = "E:/SEPSIS/RAW/SIconvert/"
BaseName ="HD214-2004_"

##
for sess in range(len(allsess)):
 
    Sessnums = allsess[sess]
    RecName = allnames[sess]
 
    files = list()
    #make file list
    for f in range(np.size(Sessnums)):
         InName = InFolder+BaseName+"%02d" % Sessnums[f]+".bin"
         list.append(files, InName)
   

   
    recording_list = list()
    for f in range(np.size(Sessnums)):
        #read axona file
        list.append(recording_list, se.AxonaRecordingExtractor(files[f]))


    ##
    MergeRec = si.concatenate_recordings(recording_list)

    ##
    TheseSigs = []
    for tet in range(np.size(Tetnums)):
       
        thistet = Tetnums[tet]
        gettet = thistet-1
        a = 4;
        chanIDs1 = gettet*a
        chanIDs1 = chanIDs1
        chanIDs2 = chanIDs1 + 4
        chanlist = range(chanIDs1,chanIDs2)
        chanlist = ["{:01d}".format(x) for x in chanlist]
       
        TheseSigs = MergeRec.get_traces(channel_ids = chanlist)
        TheseSigs = signal.decimate(TheseSigs, 4, axis = 0)
       
        thistet = str(thistet)
        outname = OutFolder+BaseName+"tt"+thistet +"_"+ RecName + ".mat"
        savemat(outname,{'Sig': TheseSigs})
       
        del TheseSigs
        TheseSigs = []
   
print('\a')