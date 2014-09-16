%% IMPORT DATA

clear all; close all;

path = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/';    % directory
date = '2014-09-16';        % date
session = '.01';
ear = '1';
cell = '13';
user = 'JS';
% include all sesion numbers (check the logdata table, columns 2 and 4)
sessions = {'132336 2';
'132358 3';
'132424 4';
'132455 5';
'132517 6';
'132543 7';
'132606 8';
'132628 9';
'132651 10';
'132733 11'}; 
statespace = '3';   % which spacefile should be used?

% Import reference data
iter_loc = 8; %logdata.data(i,iter_loc) : number of state points in the file i
alpha_loc=69;
ramp_loc1=1; %ramp_loc2=37;
spacefile = sprintf('%s%s%s%s %s%s %s%s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/State space',statespace,'.txt');
spacedata=importdata(spacefile,'\t',3);
logfile = sprintf('%s%s%s%s %s%s %s%s%s %s%s %s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',user,date,session,ear,cell,'.log');
logdata = importdata(logfile);

spacedata.data=spacedata.data; % Some spacefiles must be transposed
% Import comments
comments = logdata.textdata{2};


a = length(sessions(:,1));
for i = 1:a
file{i} = sprintf('%s%s%s%s %s%s %s%s%s %s%s %s %s %s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',user,date,session,ear,cell,'VLC',sessions{i},'.txt');
end

% Import the data 
for i = 1:a
data=importdata(file{i});
 
if i<a
    data2.data(:,:,i)=data.data;
else
    if length(data.data(:,1,1))<length(data2.data(:,1,1))
        data.data(length(data2.data(:,1,1)),:)=0;
        for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
           data.data(:,j)=zeros(1,length(data2.data(:,:)));
        end
        data2.data(:,:,i)=data.data(1:length(data2.data(:,1,1)),:);
    else
        for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
           data.data(:,j)=zeros(1,length(data2.data(:,:)));
        end
        data2.data(:,:,i)=data.data(1:length(data2.data(:,1,1)),:);
    end
end
 
end

% Import state space parameters
alpha = spacedata.data(1,:);
Fe = spacedata.data(2,:);
kv = spacedata.data(3,:);
gv = spacedata.data(4,:);
mv = spacedata.data(5,:);

Fmin = min(Fe); Fmax = max(Fe); N_F=length(unique(Fe));
kmin = min(kv); kmax = max(kv); N_k=length(unique(kv));
gmin = min(gv); gmax = max(gv); N_g=length(unique(gv));
mmin = min(mv); mmax = max(mv); N_m=length(unique(mv));

Fs = logdata.data(ramp_loc1,12);
pre = logdata.data(ramp_loc1,22)*1e-3*Fs;
pulse = logdata.data(ramp_loc1,23)*1e-3*Fs;

time = data2.data(:,1);

for i =1:a
data0(:,:,i)=data2.data(:,:,i);
end

% Import time traces
for j = 1:a
    for i = 1:(logdata.data(1,iter_loc))
        Xd(:,i,j) = data0(:,(1+i),j);
        Xo(:,i,j) = data0(:,(1+logdata.data(1,iter_loc)+i),j);
    end
end

% Extract pulses from full time traces
for j = 1:a
for i = 1:(logdata.data(1,iter_loc))
Xd_pulse(:,i,j) = Xd(1+pre:pre+pulse,i,j);
Xo_pulse(:,i,j) = Xo(1+pre:pre+pulse,i,j);
end
end

% Save extracted data
disp('Saving...');
date2= date(date~='-');
savefile = sprintf('%s%s%s%s %s%s %s%s%s%s%s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',date2,'-cell',cell,'-extracteddata.mat');
save(savefile);
disp('Finished.');
