% shortestpath() - Find shortest path to rereference Bipolar EEG
%                  of EEG/fMRI to average mastoids and data
%
% Usage:
%   >> [data_reref, A]=shortestpath(data,excludedpairs);
%
% Inputs:
%   data            - EEG dataset [channels x samples]
%   excludedpairs   - vector of bipolar pairs to exclude from rereference
%
% Outputs:
%   data_reref  - Rereferenced dataset [channels x samples]
%                   Channel assignments:
%                   1   FP1
%                   2   FP2
%                   3   AF3
%                   4   AF4
%                   5   F7
%                   6   F3
%                   7   FZ
%                   8   F4
%                   9   F8
%                   10  FC5
%                   11  FC1
%                   12  FC2
%                   13  FC6
%                   14  T7
%                   15  C3
%                   16  CZ
%                   17  C4
%                   18  T8
%                   19  CP5
%                   20  CP1
%                   21  CP2
%                   22  CP6
%                   23  P7
%                   24  P3
%                   25  PZ
%                   26  P4
%                   27  P8
%                   28  PO7
%                   29  PO3
%                   30  PO4
%                   31  PO8
%                   32  O1
%                   33  OZ
%                   34  O2
%                   35  ECG
%   A           - Rereferencing matrix, rereference data using:
%                       data_reref=A*data(1:43,:);   % Rereference data
%                       data_reref(35,:)=data(44,:); % Include ECG channel
%
% Authors: Adam Gerson (adg71@columbia.edu, 2005),
%          and Robin Goldman (rg2146@columbia.edu, 2005)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Adam Gerson and Robin Goldman
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function [data_reref, A] = shortestpath(data,excludedpairs);

if nargin<2, excludedpairs=[]; end
channels=setdiff([1:43],excludedpairs);
gain=10e3; % Amplifier gain

% Bipolar pairs
chanmap=zeros(43,36);
chanmap(1,15)=1;  chanmap(1,14)=-1; % C3-T7     1
chanmap(2,14)=1;  chanmap(2,35)=-1; % T7-LM     2
chanmap(3,35)=1;  chanmap(3,19)=-1; % LM-CP5    3
chanmap(4,19)=1;  chanmap(4,23)=-1; % CP5-P7    4
chanmap(5,23)=1;  chanmap(5,28)=-1; % P7-PO7    5
chanmap(6,28)=1;  chanmap(6,29)=-1; % PO7-PO3   6
chanmap(7,29)=1;  chanmap(7,32)=-1; % PO3-O1    7
chanmap(8,32)=1;  chanmap(8,33)=-1; % O1-Oz     8
chanmap(9,29)=1;  chanmap(9,24)=-1; % PO3-P3    9
chanmap(10,24)=1; chanmap(10,20)=-1; % P3-CP1   10
chanmap(11,25)=1; chanmap(11,20)=-1; % Pz-CP1   11
chanmap(12,20)=1; chanmap(12,15)=-1; % CP1-C3   12
chanmap(13,16)=1; chanmap(13,15)=-1; % Cz-C3    13

chanmap(14,2)=1;  chanmap(14,1)=-1; % Fp2-Fp1   14
chanmap(15,1)=1;  chanmap(15,3)=-1; % Fp1-AF3   15
chanmap(16,4)=1;  chanmap(16,2)=-1; % AF4-Fp2   16
chanmap(17,3)=1;  chanmap(17,6)=-1; % AF3-F3    17
chanmap(18,8)=1;  chanmap(18,4)=-1; % F4-AF4    18
chanmap(19,6)=1;  chanmap(19,5)=-1; % F3-F7     19
chanmap(20,9)=1;  chanmap(20,8)=-1; % F8-F4     20
chanmap(21,11)=1; chanmap(21,6)=-1; % FC1-F3    21
chanmap(22,8)=1;  chanmap(22,12)=-1; % F4-FC2   22
chanmap(23,5)=1;  chanmap(23,10)=-1; % F7-FC5   23
chanmap(24,13)=1; chanmap(24,9)=-1; % FC6-F8    24
chanmap(25,10)=1; chanmap(25,14)=-1;% FC5-T7    25
chanmap(26,18)=1; chanmap(26,13)=-1;% T8-FC6    26
chanmap(27,16)=1; chanmap(27,7)=-1; % Cz-Fz     27
chanmap(28,7)=1;  chanmap(28,11)=-1; % Fz-FC1   28
chanmap(29,12)=1; chanmap(29,7)=-1; % FC2-Fz    29
 
chanmap(30,18)=1; chanmap(30,17)=-1;% T8-C4     30
chanmap(31,36)=1; chanmap(31,18)=-1;% RM-T8     31
chanmap(32,22)=1; chanmap(32,36)=-1;% CP6-RM    32
chanmap(33,27)=1; chanmap(33,22)=-1;% P8-CP6    33
chanmap(34,31)=1; chanmap(34,27)=-1;% PO8-P8    34
chanmap(35,30)=1; chanmap(35,31)=-1;% PO4-PO8   35
chanmap(36,34)=1; chanmap(36,30)=-1;% O2-PO4    36
chanmap(37,33)=1; chanmap(37,34)=-1;% Oz-O2     37
chanmap(38,26)=1; chanmap(38,30)=-1;% P4-PO4    38
chanmap(39,21)=1; chanmap(39,26)=-1;% CP2-P4    39
chanmap(40,21)=1; chanmap(40,25)=-1;% CP2-Pz    40
chanmap(41,17)=1; chanmap(41,21)=-1;% C4-CP2    41
chanmap(42,17)=1; chanmap(42,16)=-1;% C4-Cz     42
chanmap(43,25)=1; chanmap(43,33)=-1;% Pz-Oz     43

chanmap(41,:)=-1.*chanmap(41,:); % Channel 41 is inverted

for i=1:43, pairs(i,1)=find(chanmap(i,:)==1); pairs(i,2)=find(chanmap(i,:)==-1); end

rerefmap=eye(34);
rerefmap(:,35)=-0.5; rerefmap(:,36)=-0.5; % Rereference to average mastoids
A=zeros(34,43);

% Find shortest path between electrodes
% First find distance between electrodes
eloc=readlocs('efmri36mastoids.ced','filetype','chanedit');
%eloc=readlocs('robins40.ced','filetype','chanedit');

% modified by Nelson -- August, 2006
for i=1:36
    if isempty(eloc(i).X)
        tempX(i)=0;
    else
        tempX(i)=eloc(i).X;
    end
    if isempty(eloc(i).Y)
        tempY(i)=0;
    else
        tempY(i)=eloc(i).Y;
    end
    if isempty(eloc(i).Z)
        tempZ(i)=0;
    else
        tempZ(i)=eloc(i).Z;
    end
end
xmat=repmat([tempX],[36 1]);
xdist=xmat-xmat';
ymat=repmat([tempY],[36 1]);
ydist=ymat-ymat';
zmat=repmat([tempZ],[36 1]);
zdist=zmat-zmat';
dist=sqrt(xdist.^2+ydist.^2+zdist.^2); % Distance between electrodes


% Construct distance matrix for bipolar electrode graph
% First set d(i,j) to 1 if arc exists between electrodes i and j
% and set d(i,j) to Inf if arc does not exist
d(1,:,:)=eye(36);
for i=1:36,
   d(1,i,pairs(find(pairs(:,1)==i),2))=1;
   d(1,i,pairs(find(pairs(:,2)==i),1))=1;
end
d(find(d==0))=Inf;

% Remove excluded pairs
for i=1:length(excludedpairs),
   d(1,pairs(excludedpairs(i),1),pairs(excludedpairs(i),2))=Inf; %%problem here
   d(1,pairs(excludedpairs(i),2),pairs(excludedpairs(i),1))=Inf;
end

% Now set d(i,j) to dist(i,j) if arc exists between electrodes i and j
d(1,:,:)=squeeze(d(1,:,:)).*dist;

% Construct initial predecessor matrix
p(1,:,:)=repmat([1:36]',[1 36]);
for i=1:36,
   p(1,i,i)=NaN;
end

for k=2:37,
   for i=1:36,
       for j=1:36,
           d(k,i,j)=min(d(k-1,i,j),d(k-1,i,k-1)+d(k-1,k-1,j));
           if d(k,i,j)~=d(k-1,i,j),
               p(k,i,j)=p(k-1,k-1,j);
           else
               p(k,i,j)=p(k-1,i,j);
           end
       end
   end
end

% Final solution for shortest paths between all pairs of electrodes
p=squeeze(p(37,:,:)); d=squeeze(d(37,:,:));

% Find rereferencing sequences of bipolar pairs
% Since there are three possible paths between 3 nodes,
% first search for shortest path that allows electrode to be rereferenced
% to average mastoids.  If these three shortest paths do not have a
% rereferencing solution, then use the general solution.
% Ideally, rather than using the general solution it would be better to try
% the next shortest path.
% The shortest path actually appears to consistently yield a
% rereferencing solution.  However see kthshortestpath.m if the need ever
% arises for this solution.  ADG

for channel=1:34,
   pselect=1; lastwarn(''); lastid='MATLAB:rankDeficientMatrix';
   while strcmp(lastid,'MATLAB:rankDeficientMatrix')&(pselect<5)
       if pselect<4,
           pairi=findpath(channel,pselect,d,p,pairs);
           A(channel,:)=0;
           A(channel,pairi)=rerefmap(channel,:)/chanmap(pairi,:);
       end
       if pselect==4,
           A(channel,:)=0;
           A(channel,channels)=rerefmap(channel,:)/chanmap(channels,:); lastwarn('');
       end
       [lastmsg lastid]=lastwarn; lastwarn('');
       pselect=pselect+1;
       %if strcmp(lastid,'MATLAB:rankDeficientMatrix'),
       %    fprintf(['Channel ' eloc(channel).labels ' rank deficient, recalculating\n']);
       %    A(channel,channels)=rerefmap(channel,:)/chanmap(channels,:); lastwarn('');
       %end
   end
%    fprintf(['Channel ' eloc(channel).labels ', path = %d\n'],pselect-1);
end % channel


% Round to 6 decimal places
A=round(A.*10^6)./(10^6); % This is the rereferencing matrix

data_reref=A*data(1:43,:);   % Rereference data
data_reref(35,:)=data(43,:); % Include ECG channel
%data_reref(34,:)=data(43,:); % changed by ian 7/29/13
% commented out by Jen because preprocessing code adjusts for gain:
%data_reref=1e6.*data_reref./gain; % Convert to microvolts 

function pairindx=findpath(channel,pselect,d,p,pairs)
pathd=[d(channel,35) d(channel,36) d(35,36)];
paths=[channel 35; channel 36; 35 36];
[sortp,sorti]=sort(pathd);
if pselect==1, sorti=sorti([1 2]); end
if pselect==2, sorti=sorti([1 3]); end
if pselect==3, sorti=sorti([2 3]); end
pairindx=[];
tmp=pairsequence(p,pairs,paths(sorti(1),1),paths(sorti(1),2)); % LM -> channel
pairindx=[pairindx; tmp];
tmp=pairsequence(p,pairs,paths(sorti(2),1),paths(sorti(2),2)); % RM -> channel
pairindx=[pairindx; tmp];
pairindx=unique(pairindx);

function [pairindx]=pairsequence(p,pairs,node1,node2)
pairindx=[];
pindx=node1;
i=1;
while pindx(i)~=node2,
   i=i+1;
   pindx(i)=p(node2,pindx(i-1));
   tmp=find(pairs(:,1)==pindx(i));
   if ~isempty(find(pairs(tmp,2)==pindx(i-1))),
       pairindx(i-1,1)=tmp(find(pairs(tmp,2)==pindx(i-1)));
   else
       tmp=find(pairs(:,2)==pindx(i));
       pairindx(i-1,1)=tmp(find(pairs(tmp,1)==pindx(i-1)));
   end
end

