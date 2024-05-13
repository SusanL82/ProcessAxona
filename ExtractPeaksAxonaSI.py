
import spikeinterface.extractors as se
import spikeinterface as si
import numpy as np
from probeinterface import read_prb
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.preprocessing import bandpass_filter
from scipy.io import savemat
from tqdm import tqdm


InFolder = "E:/AxonaToNLXtest/"
OutFolder = "E:/AxonaToNLXtest/pos/"
Probepath ="C:/Users/leemburg/Desktop/OEphystest/"
ChanList = 'Tetlist.txt' # text file listing good and bad channels
BaseName = 'HD263-2811_'
TetList = [2,4,6] #analyse only these tetrodes (1-based, as on drive)
numsess = 3 #number of sessions to read (e.g: numsess = 4 reads bin 01 to 04)

spike_thresh = 5 # detection threshold for spike detection is spike_thresh* signal SD
spike_sign = 'pos' #detect peaks. can be: 'neg', 'pos', 'both'

# !!! need to manually edit the .set files so that all channels get read. 

##############################################################################

tetgrouping = np.array([0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,
                        9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,
                        15,15,15,15,15])

# read bad channel list, set bad channels to 17
tetlist = np.loadtxt(InFolder + '/' + ChanList,
                 delimiter=",")

for tetnum in range(16):
        thistet = np.where(tetgrouping==tetnum)
        thistet = np.array(thistet)
        thesewires = tetlist[tetnum,1:5]
        badwires = np.where(thesewires==0)
        badwires = np.array(badwires)
        badchans = thistet[0][badwires]
        tetgrouping[badchans]=17



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

all_chan_ids = MergeRec_f.get_channel_ids()

# select a tetrode
for tetnum in TetList:

    tet_chan_ids = all_chan_ids[np.where(tetgrouping == tetnum-1)]
    tetname = TetList[tetnum-1]

    if np.size(tet_chan_ids)>2:
        
        # 4-wire tetrode
        if np.size(tet_chan_ids) == 4:
           new_chans = [0,1,2,3]
           probename = "tet4_probe.prb"
        
        # 3-wire tetrode
        if np.size(tet_chan_ids) == 3:
            new_chans = [0,1,2]
            probename = "tet3_probe.prb"
        
    myprobe = read_prb(Probepath + "/" + probename)

    #select channels and add probe
    thistet = MergeRec_f.channel_slice(tet_chan_ids, renamed_channel_ids=new_chans)
    thistet = thistet.set_probegroup(myprobe)

    # preprocess (filter)
    thistet_f = bandpass_filter(thistet, freq_min=600, freq_max=6000)

    #detect peaks
    detectradius = 30 #tetrode channel group is 10x10um square, this radius should capture all spikes in this channelgroup
    #2 versions of function from 2 versions of spikeinterface. Yuck.
    #TetPeaks = detect_peaks(thistet_f, method='locally_exclusive', peak_sign=spike_sign, detect_threshold=spike_thresh, exclude_sweep_ms=0.1, local_radius_um=detectradius, noise_levels=None,)
    TetPeaks = detect_peaks(thistet_f, method='locally_exclusive', peak_sign=spike_sign, detect_threshold=spike_thresh, exclude_sweep_ms=0.1, radius_um=detectradius, noise_levels=None,)


    # get waveforms. must be 32 samples long, I'm setting peak at sample 8 (same as NLX)
    prepeak = 8
    postpeak = 24

    print('extracting waveforms')

    allWFs = np.zeros([32, 4, len(TetPeaks['sample_index'])], dtype='int16')
    for p in tqdm(range(len(TetPeaks['sample_index'])), desc="collecting waveforms"):
        sf = TetPeaks['sample_index'][p] - prepeak
        ef = TetPeaks['sample_index'][p] + postpeak

        thisWF = thistet_f.get_traces(segment_index=None, start_frame=sf, end_frame=ef)

        #write only complete spike waveforms (might skip last spike if too short)
        if np.size(thisWF, 0) == 32:
           allWFs[:, :, p] = thisWF

        del thisWF

   #save peaks to mat file (for later processing with Mat2NlxSpike to generate .ntt file)
    print('saving to mat file')
    outname = OutFolder+"tt"+str(tetname) + '.mat'
    savemat(outname, {"Timestamps":TetPeaks['sample_index'], "Spikes": allWFs})        
        
