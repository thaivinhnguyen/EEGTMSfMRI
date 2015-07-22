%% Load in data set
EEG=pop_loadset();
channelsToKeep=[];
EEG.data(setdiff(1:34,channelsToKeep),:)=0;
%% Run ICA
EEG=pop_runica(EEG,'fastica');
%% Visualization using EEGLAB
ALLEEG=EEG;
eeglab redraw
%% Choose ICA Component for Frontal Alpha
icaComp=3;
%% Extract Timeseries
ts=EEG.icaact(icaComp,:);
%% Choose upper and lower filter bands
lbound=8;
hbound=12;
%% Filter Around alpha
tsFilt=eeg_filt(double(ts),256,[lbound hbound]);
%% plot time series
times=linspace(EEG.xmin,EEG.xmax,EEG.pnts);
figure;plot(times,ts,'r',times,tsFilt,'b')
%% run hilbert transform on filtered data
[X]=hilbert(tsFilt);

%% plot amplitude
hold on
plot(times,abs(X),'g','LineWidth',3)
%% caculate phase
phi=angle(X);
%% plot phase
plot(times,phi,'k','LineWidth',3)

%% Collect fake TMS pulses
phis=phi(1500:1500:end);

%% Bin data
indx_neg_pi_pi2=find(phis>=-pi & phis<-pi/2);
indx_neg_pi2_0=find(phis>=-pi/2 & phis<0);
indx_pos_0_pi2=find(phis>=0 & phis<pi/2);
indx_pos_pi2_pi=find(phis>=pi/2 & phis<pi);
data=[mean(phis(indx_neg_pi_pi2)),mean(phis(indx_neg_pi2_0)),mean(phis(indx_pos_0_pi2)),mean(phis(indx_pos_pi2_pi))];
figure;bar(data)
