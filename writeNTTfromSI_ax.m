%% export spikeinterface detected peaks to an .NTT file
% uses Mat2NlxSpike v 6.0.0

%% user inputs
% Inpath = 'D:\SEPSIS\ntt_test'; % path with input file
% Outpath = 'D:\SEPSIS\ntt_test'; %path for storing the generated .ntt
% Filename = 'HD190-1511_tt2.mat'; %name of input file, this is the .mat generated
% by spikeinterface with ExtractPeaksAxonaSI,py

function [Filename] = writeNTTfromSI_ax(Inpath, Outpath, Filename)

%% load SI data
%get waveforms and timestamps
load([Inpath,'\',Filename])

%% write .ntt

Outname = strsplit(Filename,'.mat');
Outname = Outname{1};
NTTname = [Outpath,'\',Outname,'.ntt'];

% FieldSelectionFlags(1): Timestamps (1xN vector of timestamps, ascending
% order
% FieldSelectionFlags(2): Spike Channel Numbers
% FieldSelectionFlags(3): Cell Numbers (here: 0, no cells sorted yet)
% FieldSelectionFlags(4): Spike Features (8xN integer vector of features
% from cheetah: peaks for 4 channels and valley for 4 channels.
% FieldSelectionFlags(5): Samples 32x4xN integer matrix with the datapoints
% (waveform) for each spike for all 4 channels.
% FieldSelectionFlags(6): Header

AppendToFileFlag = 0; %new file will be created or old file will be overwritten
ExportMode = 1; %export all
FieldSelectionFlags = [1,1,1,1,1,0];

numspikes = numel(Timestamps);
ScNumbers = zeros(1,numspikes); %set to 0 (cheetah also does that)
CellNumbers = ScNumbers; %all cells in cluster 0

Features = nan(8,numel(Timestamps));
for s = 1:numel(Timestamps)
    Features(1:4,s) = max(Spikes(:,:,s));
    Features(5:8,s) = min(Spikes(:,:,s));
end

disp(['exporting to ', Outpath])
Mat2NlxSpike(NTTname, AppendToFileFlag, ExportMode, [], FieldSelectionFlags, Timestamps, ScNumbers, CellNumbers, Features, Spikes)
disp(['created ',Filename,'.ntt'])

end