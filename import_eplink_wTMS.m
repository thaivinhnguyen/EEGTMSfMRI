function EEG=import_eplink_wTMS(filename,cap)


%% read from raw binary file
[EEG_data,EventChannel,~]=read_amp_binaryv3(filename,cap);

%% removal gradient artitifact
[EEG_GA_Removed]=gradient_removal(EEG_data,2,1);


%% Check Event Channel is Same Size as data chans
if size(EEG_GA_Removed,2)~=size(EventChannel,2)
    EventChannel=EventChannel(1:size(EEG_GA_Removed,2));
end

%% Add in fake events for TMS
EventChannel(1:5000:end)=100;

%% Combine EventChannel and Data Channels
EEG_GA_Total=[EEG_GA_Removed;EventChannel'];

% save([fullfile(pathnames{1},filename) '.mat'],'EEG_data');
%
% save([fullfile(pathnames{2},filename) '_Gradient_Removed.mat'],'EEG_GA_Total');
clear EEG_GA_Removed
[a b]=fileparts(filename);

ALLEEG=[];
setname = [b '_Gradient_Removed'];
fname = setname;
fprintf(['\n--------------------------------------------------------------' ...
    '\n   Importing %s to EEGLAB' ...
    '\n--------------------------------------------------------------\n'],fname)
EEG=pop_importdata('dataformat','array','data',EEG_GA_Total,'srate',488,'setname', setname);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

EEG = pop_chanevent(EEG,EEG.nbchan, 'edge', 'leading', ...
    'edgelen', 1, 'delchan', 'on', 'delevent', 'on', 'duration','off', 'nbtype', NaN);
EEG = eeg_checkset(EEG);

EEG = pop_select(EEG, 'nochannel', 44+1:EEG.nbchan);
EEG = eeg_checkset(EEG);
