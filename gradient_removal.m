function [EEG_GA_Removed]=gradient_removal(EEG_data,tr,medFilt)

%%
% 1. Epoch by the master 000/ each epoch should start with 3 zeros/by
% channel
% 2. Average across epoch within channel
% 3. Subjract the average within channel
% 4. see whats remaining

%% NEED TO DO
% Unhardcode the TR length in datapoints
% figure out if we should use filter or filtfilt
% Found bug in one resting state scan. First 4 datapoints are 0

%%
if nargin<3
    medFilt=1;
end
fs=488;
EEG_GA_Removed=zeros(size(EEG_data));
% [hpnum,hpdenom]=butter(2,1/fs*2,'high');
shiftSamples=tr*fs-2;
fprintf('Removing Gradient from Channel ');
for chan=1:size(EEG_data,1)
    fprintf('%d ',chan)
    EEG_data_tmp=detrend(EEG_data(chan,:));
    if sum(EEG_data(chan,:))~=0
        chan_zeros=find(EEG_data(chan,:)==0);
        %check if the first 4 samples are 0;
        sumNum=sum(chan_zeros(1:4));
        if sumNum(1)==10
            chan_zeros(1)=[];
        end
        count=1;
        first_zero=[];
        for z=1:length(chan_zeros)-2,
            if sum(EEG_data(chan,chan_zeros(z):chan_zeros(z)+2))==0
                first_zero(count)=chan_zeros(z);count=count+1;
            end
        end
        for TR=1:length(first_zero)
            tmp(TR,:)=squeeze(EEG_data_tmp(first_zero(TR)-1:first_zero(TR)+shiftSamples));
        end
%         chan_GA=bsxfun(@minus,squeeze(tmp(chan,:,:)),mean(squeeze(tmp(chan,:,:)),1));
% chan_GA=reshape(chan_GA,1,prod(size(chan_GA)));
        no_GA_mean=mean(squeeze(tmp(:,:)),1);
        
        for TR=1:length(first_zero)
            EEG_GA_Removed(chan,first_zero(TR)-1:first_zero(TR)+shiftSamples)=...
                EEG_data_tmp(first_zero(TR)-1:first_zero(TR)+shiftSamples)-no_GA_mean;
        end
        if medFilt
        EEG_GA_Removed(chan,:)=stuff(EEG_GA_Removed(chan,:));
        end
    end
    
    clear tmp
end

fprintf('\nDone\n')
% EEG_GA_Removed(end,:)=EEG_data(end,:);