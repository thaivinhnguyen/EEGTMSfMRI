% [data_noBCG, BCGartifact] = removeBCG(data,fs,eegchan,includedsamples,BCG_npc);
%
% Remove BCG from EEG/fMRI data using max power method
%
% Nelson Wei (cw2089@columbia.edu) on 11/18/2006
% Jennifer Walz (jw2552@columbia.edu) extracted just the BCG part of code
% on 2/1/10
% number of BCG components is now an input to the function - JMW 7/13/10

% # Channels: 43 EEG, 2 EOG, 2 ECG, 2 Events;  # Electrodes: 34 EEG
% Sampling Freq: 1 kHz



function [data, bcgartifact, tmp] = removeBCG(data,fs,eegchan,includesamples,BCG_npc)

if isempty(eegchan), eegchan = 1:size(data,1); end
if isempty(includesamples), includesamples = 1:size(data,2); end
if isempty(BCG_npc), BCG_npc = 2; end

% LPF @ 4Hz to find BCG
data_lp = eegfilt(data(eegchan,includesamples),fs,0,4);

%BCG_npc = 2; % number of principle components 

% PCA (by finding the eigenvectors of the covariance matrix)
cov_data = data_lp*data_lp';  % covariance of the lowpassed EEG data
[vb,tmp] = eig(cov_data);     % find eigenvectors

% JEN ADD:
% sort the eigenvalues just in case Matlab isn't doing it right
[tmp,idx] = sort(diag(tmp),'ascend'); 
% sort the eigenvectors in ascending order 
vb = vb(:,idx);
% fprintf('\nThe sorted eigenvalues are:\n')
% tmp
% fprintf('\nThe sorted eigenvectors are:\n')
% vb
% fprintf('\nThe eigenvectors were sorted using these indices:\n')
% idx

normmaxcomp=norm(vb(:,end)); % find norm of first principal component

% divide each principle component by length of the 1st component
for i=1:BCG_npc, 
    Abcg(:,i)=vb(:,end-(i-1))./normmaxcomp; 
end

% use least squares to solve?
Vbcg=inv(Abcg'*Abcg)*Abcg';
bcgartifact(eegchan,includesamples) = Abcg*Vbcg*data(eegchan,includesamples);
% subtract out the artifact
data(eegchan,includesamples) = (eye(length(eegchan))-Abcg*Vbcg)*data(eegchan,includesamples);



