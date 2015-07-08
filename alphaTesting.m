%%
close all
clear all
EEG=pop_loadset();
[spectra,freqs,speccomp,contrib,specstd] = ...
    spectopo(EEG.data, 0, EEG.srate,'freqrange',[2 25]);
nonAlphaFreqs=(freqs>=5 & freqs<=7.5) | (freqs>=12.5 & freqs<=15);
AlphaFreqs=freqs(freqs>=7.5 & freqs <=12.5);
[alphaPow,I]=max(spectra(:,freqs>=7.5 & freqs <=12.5),[],2);
nonPow=mean(spectra(:,nonAlphaFreqs),2);





I=10.^(alphaPow/10)./10.^(nonPow/10);
Ipow=10*log(I);
figure;subplot(1,2,1)
topoplot(alphaPow,EEG.chanlocs),caxis([min(alphaPow) max(alphaPow)])

subplot(1,2,2)
topoplot(Ipow,EEG.chanlocs),caxis([min(Ipow) max(Ipow)])
%  [smoothdataAlpha] = eegfilt(EEG.data,EEG.srate,7,13,0,[],0);
%  [smoothdataNonAlpha]= eegfilt(EEG.data,EEG.srate,7,13,0,[],1);

t=5;%5 second chunks
tnpts=5*EEG.srate;
[chan,npts]=size(EEG.data);
chunks=round(npts/tnpts)-2;

for j=1:chunks
    
    [spectra,freqs,speccomp,contrib,specstd] = ...
        spectopo(EEG.data(:,tnpts*(j)+1:(j+1)*tnpts),0, EEG.srate,'freqrange',[2 25],'plot','off');
    nonAlphaFreqs=(freqs>=5 & freqs<=7.5) | (freqs>=12.5 & freqs<=15);
    AlphaFreqs=freqs(freqs>=7.5 & freqs <=12.5);
    [alphaPow2(:,j),I]=max(spectra(:,freqs>=7.5 & freqs <=12.5),[],2);
    freqDist(:,j)=AlphaFreqs(I);
    nonPow2(:,j)=mean(spectra(:,nonAlphaFreqs),2);
end

I2=10.^(alphaPow2./10)./10.^(nonPow2./10);
I2pow=10*log(I2);

figure;subplot(1,3,1)
topoplot(mean(alphaPow2,2),EEG.chanlocs),caxis([min(mean(alphaPow2,2)) max(mean(alphaPow2,2))])

subplot(1,3,2)
topoplot(mean(I2pow,2),EEG.chanlocs),caxis([min(mean(I2pow,2)) max(mean(I2pow,2))])

subplot(1,3,3)
topoplot(std(I2pow,[],2),EEG.chanlocs),caxis([min(std(I2pow,[],2)) max(std(I2pow,[],2))])

