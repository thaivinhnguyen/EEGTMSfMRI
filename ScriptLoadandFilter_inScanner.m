%% Script to load and preprocess EEG-fMRI data
% add script location to path
P = mfilename('fullpath');
[direc,~]=fileparts(P);
addpath(direc)

%filter settings
bandpass_bounds = [0.5 100];
notch_bounds = [59 61];
% Set filter bounds
highpass_cutoff = bandpass_bounds(1);
lowpass_cutoff = bandpass_bounds(2);

% make sure eeglab is running
eeglab_is_running = evalin('base','exist(''ALLEEG'',''var'')'); % if ALLEEG is a variable in the base class, eeglab is running.
if ~eeglab_is_running
    eeglab;
end


%% Need to change to where the loc file directory
% electrodeloc_dir = '/Volumes/lorentz2/BASEBALL/';
pathnames=pwd; %save files to current directory
prefix='restingstate';
eegSession=[];

filename_base = sprintf('%s%d',prefix,eegSession);
final_filename = filename_base;

raw_fs = 488; % sampling frequency of raw data
new_fs = 256;
lpf_useIIR = 0;
lpf_order = [];
eventchans=47;
ref=1:1:46;
cap=1;

EEG = import_eplink(final_filename,cap);

EEG = eeg_checkset( EEG );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',[filename_base '-raw'],'gui','off');


%% Low-Pass Filter the data
if ~isnan(lowpass_cutoff) && highpass_cutoff < lowpass_cutoff
    [EEG, ~, lpf] = pop_eegfilt_returnfilter( EEG, 0, lowpass_cutoff, lpf_order, 0, lpf_useIIR); %low-pass, 100Hz cutoff
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-lowpass'],'overwrite','on','gui','off');
    EEG = eeg_checkset( EEG );
end

%% Resample data to new_fs Hz
if new_fs < raw_fs
    EEG = pop_resample_myfilter(EEG,new_fs,0); % do not apply filter
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-resampled'],'overwrite','on','gui','off');
    EEG = eeg_checkset( EEG );
end

%% Filter the data further
if ~isnan(highpass_cutoff) && highpass_cutoff < lowpass_cutoff
    EEG = pop_eegfilt( EEG, highpass_cutoff, 0, [], [0]); %high-pass, 0.5Hz cutoff
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-highpass'],'overwrite','on','gui','off');
    EEG = eeg_checkset( EEG );
end
if notch_bounds(2)>notch_bounds(1)
    EEG = pop_eegfilt( EEG, notch_bounds(1), notch_bounds(2), [], [1], [1]); %notch, 59-61Hz, FFT
else
    warning('DJ:NotchBounds','Skipping Notch Filter!')
end
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename_base '-filtered'],'overwrite','on','gui','off');

%% Save the data
EEG = pop_saveset( EEG, 'filename',final_filename, 'filepath',pathnames);

%% BCG Removal
[data,bcgartifact]=removeBCG(EEG.data,EEG.srate,[],[],[]);
EEG.data=data;

EEG = pop_saveset(EEG,'filename', [final_filename '_noBCG'],'filepath',pathnames);

%% Reference Data
% EEG = pop_loadset([final_filename '.set']);
%% re-reference
if ~exist(fullfile(pathnames,[final_filename,'_eps.mat']),'file')
    pop_eegplot(EEG)
    %Remove bad channels
    eps=str2num(input('Excluded Pairs:','s'));
    
    [data_reref,A]=shortestpath(EEG.data,eps);
    data_reref=data_reref(1:34,:);
    save(fullfile(pathnames,[final_filename,'_eps']),'eps','A');
else
    load(fullfile(pathnames,[final_filename,'_eps']));
    data_reref=A*EEG.data(1:43,:);
end
EEG.data=data_reref;
EEG=eeg_checkset(EEG);
EEG.chanlocs=readlocs('efmri34.ced');
EEG=eeg_checkset(EEG);
EEG=pop_saveset(EEG,'filename',[EEG.filename(1:end-4) '_noEPs.set'],'filepath',pathnames);


