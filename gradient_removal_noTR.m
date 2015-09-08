function [EEG_GA_Removed]=gradient_removal_noTR(EEG_data,numFilt,medFilt)
% EEG_data-- a matrix of channels x timepoints
% numFilt-- number of TRs to use in moving average (if 0 uses all TRs)
% medFilt-- 1 to apply median filter

%%
% 1. Epoch by the master 000/ each epoch should start with 3 zeros/by
% channel
% 2. Average across epoch within channel
% 3. Subjract the average within channel
% 4. see whats remaining
% 5. Unhardcode the TR length in datapoints--Fixed
% 6. Using regress instead of subtraction (doesn't change much but handles
% drifts in data better
%% NEED TO DO
% Found bug in one resting state scan. First 4 datapoints are 0

%%
if nargin<2
    numFilt=0;
end
if nargin<3
    medFilt=1;
end
EEG_GA_Removed=zeros(size(EEG_data));
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
        shiftSamples=round(mean(diff(first_zero)))-2;
        tmp=zeros(length(first_zero),shiftSamples+2);
        for TR=1:length(first_zero)
            tmp(TR,:)=squeeze(EEG_data_tmp(first_zero(TR)-1:first_zero(TR)+shiftSamples));
        end
        %         chan_GA=bsxfun(@minus,squeeze(tmp(chan,:,:)),mean(squeeze(tmp(chan,:,:)),1));
        % chan_GA=reshape(chan_GA,1,prod(size(chan_GA)));
        if numFilt<5 || isempty(numFilt)
            no_GA_mean=mean(squeeze(tmp(:,:)),1);
            
            for TR=1:length(first_zero)
                [~,~,R]=regress(EEG_data_tmp(first_zero(TR)-1:first_zero(TR)+shiftSamples)',no_GA_mean');
                EEG_GA_Removed(chan,first_zero(TR)-1:first_zero(TR)+shiftSamples)=R';
            end
        else
            
            
            for TR=1:length(first_zero)
                if TR<=(numFilt/2)
                    trs=setdiff([1:numFilt],TR);
                elseif TR+numFilt/2>length(first_zero) 
                    trs=setdiff([length(first_zero)-numFilt:length(first_zero)],TR);
                else
                    trs=setdiff(TR-numFilt/2:TR+numFilt/2,TR);
                end
                no_GA_mean=mean(tmp(trs,:),1);
                % Get the residuals
                [~,~,R]=regress(EEG_data_tmp(first_zero(TR)-1:first_zero(TR)+shiftSamples)',no_GA_mean');
                EEG_GA_Removed(chan,first_zero(TR)-1:first_zero(TR)+shiftSamples)=R';
            end
            
        end
        
        
        if medFilt
            EEG_GA_Removed(chan,:)=stuff(EEG_GA_Removed(chan,:));
        end
    end
    
    clear tmp
end

fprintf('\nDone\n')
