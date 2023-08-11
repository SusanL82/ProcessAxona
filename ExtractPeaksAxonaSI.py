# -*- coding: utf-8 -*-

import spikeinterface as si
import spikeinterface.extractors as se
from scipy.io import savemat
import numpy as np
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.preprocessing import bandpass_filter
from tqdm import tqdm
from probeinterface import Probe

InFolder = "D:/SEPSIS/"
OutFolder = "D:/SEPSIS/ntt_test/"
BaseName ="HD190-1511_"
ChanList = 'tetlist.txt' # text file listing good and bad channels, must be in InFolder

numsess = 2 #number of sessions/ .bin files
Tetnums = [2,4,5,6,7,8,9,11] #tetrodes to process (1-based, same as on the rat/rec)
spike_thresh = 5 # detection threshold for spike detection is spike_thresh* signal SD
spike_sign = 'neg' #detect negative peaks only. can be: 'neg', 'pos', 'both'
spike_sweep = 0.1 #exclude spikes detected within x ms. Cheetah uses 0.75ms
spike_radius = 30 #local neighborhood for exclusive spike detection (in this case, the whole tetrode group where tetrodes are 10um apart)


# !!! need to manually edit the .set files so that all channels get read. 

##############################################################################

# get bad channel list
tetlist = np.loadtxt(InFolder + '/' + ChanList,
                 delimiter=",")

# make file list
allsess = range(numsess) #count sessions zero-based
files = list()
for f in allsess:
   thissess = allsess[f]+1 #1-based session number
   InName = InFolder+BaseName+"%02d" % thissess+".bin"
   list.append(files, InName)


# make list of SI extractor recordings
recording_list = list()
for f in allsess:
   list.append(recording_list, se.AxonaRecordingExtractor(files[f]))

# merge recording objects 
MergeRec = si.concatenate_recordings(recording_list)

# add highpass filter
MergeRec_f = bandpass_filter(MergeRec, freq_min=300.0, freq_max=6000.0, margin_ms=5.0, dtype=None)


# detect peaks per tetrode and get trace segments for waveforms
for tet in range(np.size(Tetnums)):
       
        thistet = Tetnums[tet]
        tetchans = tetlist[thistet-1,1:5]   
        
        print('processing TT'+str(thistet))
        
        # get channel IDs for good channels
        gettet = thistet-1
        a = 4;
        chanIDs1 = gettet*a
        chanIDs2 = chanIDs1 + 4
        chanlist = range(chanIDs1,chanIDs2)
        chanlist = chanlist * tetchans
        
        goodchans = np.where(tetchans==1)
        goodchans = np.array(goodchans)
        chanlist = chanlist[goodchans]
        
        #make tetrode layout
        xcoords = [0,0,10,10]
        ycoords = [0,10,0,10]
        #make probe for these channels
        positions = np.zeros((np.size(chanlist), 2))
        for p in range(np.size(chanlist)):
            x = xcoords[p]
            y = ycoords[p]
            positions[p] = x, y
        probe = Probe(ndim=2, si_units='um')
        probe.set_contacts(positions=positions, shapes='circle', shape_params={'radius': 5})
        
        
        chanlist = ["{:01d}".format(x) for x in chanlist]
       
        TheseSigs = MergeRec_f.channel_slice(channel_ids = chanlist)
        TheseSigs = TheseSigs.set_probe(probe)
        TetPeaks = detect_peaks(TheseSigs, method='locally_exclusive', pipeline_nodes = None, gather_mode='memory', folder=None, names=None, 
                                peak_sign=spike_sign, detect_threshold= spike_thresh, exclude_sweep_ms = spike_sweep, radius_um = spike_radius)
                 #TetPeaks has fields: sample_index, channel_index, amplitude, segment_index
        
        
        # get waveforms. must be 32 samples long, I'm setting peak at sample 8 (same as NLX)
        prepeak = 8 
        postpeak = 24
       
        print('extracting waveforms')
        
        allWFs = np.zeros([32,4,len(TetPeaks['sample_index'])], dtype = 'int16') 
        for p in tqdm(range(len(TetPeaks['sample_index'])), desc="collecting waveforms"):
            sf = TetPeaks['sample_index'][p] - prepeak
            ef = TetPeaks['sample_index'][p] + postpeak
        
            thisWF = TheseSigs.get_traces(segment_index = None, start_frame = sf, end_frame = ef)
           
            #write only complete spike waveforms (might skip last spike if too short)
            if np.size(thisWF,0) == 32:
                allWFs[:,:,p] = thisWF           
            
            del thisWF
       
        #save peaks to mat file (for later processing with Mat2NlxSpike to generate .ntt file)
        print('saving to mat file')
        outname = OutFolder+BaseName+"tt"+str(thistet) +'.mat'
        savemat(outname,{"Timestamps":TetPeaks['sample_index'], "Spikes": allWFs})e_index'], "Spikes": allWFs})
        
        
