function [EEG_data,EventChannel,AmpInfo]=read_amp_binaryv3(filename,CapType)
%%
% Version 1 of new amp channel data read
% first 64 bits are amplifier info, then 128 bit24 data points are read in
% and then one byte is read for the event channel and repeats
%Captype indicates what channels we use for the currect cap
if nargin<2
    [filename, pathname] = uigetfile;
end

if nargin<3
    CapType=1;
end

if CapType==1
     eeg_channels=[1:48];
     bo_box=[127,128];
%     channels=[eeg_channels,bo_box];
    % Correct Channel Placements
    % OLD CHANNEL PLACEMENTS AND WRONG 
    channels=[33:45 1:31 46 127 128]; 
%     channels=[33:45
else
    channels=[1:1:128];
    bo_box=[];
end
% Read the contents back into an array
fid = fopen(fullfile(filename));

[AmpInfo]=fread(fid,64,'ubit1','ieee-be');
AmpInfo=reshape(AmpInfo,8,8);


%Check Gains are the Same within eeg_channels
if CapType==1
    %we know first 6 cards are on, so check all 6 are equal
    if sum(sum(AmpInfo([1,5],1:3)))==6
        gains=sum([AmpInfo(3:4,1:3)';AmpInfo(7:8,1:3)']);
        
        %check that you have either a sum of 6 or 0;
        if sum(gains)~=0 && sum(gains)~=6
            error('AmpInfo Wrong, Not All cards have same gain')
        else
            if gains(1)==0 && gains(2)==0
                gain=400;
            elseif gains(1)==0 && gains(2)==6
                gain=10E3;
            elseif gains(1)==6 && gains(2)==0
                gain=2E3;
            else
                error('AmpInfo Wrong')
                
            end
        end
        
        
    else
        error('AmpInfo Wrong,Not all cards are on for this Cap!');
    end
    
    
end

if ~isempty(bo_box)
    %check if card is on
    if AmpInfo(1,8)==1
        if AmpInfo(3,8)==0 && AmpInfo(4,8)==0
            boGain=400;
        elseif AmpInfo(3,8)==1 && AmpInfo(4,8)==0
            boGain=2E3;
        elseif AmpInfo(3,8)==0 && AmpInfo(4,8)==1
            boGain=10E3;
        end
    else
        error('Card 16 was not Turned on')
    end
    
end

%ADC Calculations
Vref=2.5;
bitReso=2^23-1; %24 bit
microV=1E6;


fprintf('Reading in Channel ');
%%
for k=1:128
    fprintf('%d ',k)
    
    tmp=fread(fid,inf,'bit24',127*24+8,'ieee-be');
    %check if tmp has same amount of data, might not be the same if daq was
    %stopped mid recording. If different, pad ending with 0s
    
    
    
    if k==1
        EEG_data=zeros(128,numel(tmp));
    else
        if numel(tmp)~=size(EEG_data,1)
            tmp=[tmp;zeros(size(EEG_data,2)-numel(tmp),1)];
        end
    end
    if CapType
        if sum(ismember(eeg_channels,k))
            tmp=tmp*Vref/gain/bitReso*microV; %Scale EEG channels by correct gain
        elseif sum(ismember(bo_box,k))
            tmp=tmp*Vref/boGain/bitReso*microV; %Scale boBox by correct gain
        end
    else
        
    end
    
    EEG_data(k,:) = tmp';
    fseek(fid,8+k*3,'bof');
    % EventChannel(k)=tmp2;
    
    
end

%get event channel
fseek(fid,8+(128*3),'bof');
EventChannel=fread(fid,inf,'uint8',128*3,'ieee-be');


fclose(fid);

fprintf('\nDone\n');
EEG_data=EEG_data(channels,:);

