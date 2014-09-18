%% IMPORT DATA

clear all; close all;

path = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/';    % directory
date = '2014-09-17';        % date
session = '.01';
ear = '1';
cell = '4';
user = 'JS';
% include all sesion numbers (check the logdata table, columns 2 and 4)
sessions = {'133007 11';
'133036 12'}; 
statespace = '6';   % which spacefile should be used?

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
savefile = sprintf('%s%s%s%s %s%s %s%s%s%s%s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',date2,'-cell',cell,'-extracteddata-C.mat');
save(savefile);
disp('Finished.');


%% Calculate FFT / PSD for each pulse

clear Xpsd Xsinepsd fsinepsd fpsd Xmtm fmtm

% Center signals
Nsm = 10;
for j = 1:a
    for i = 1:(logdata.data(1,iter_loc))
        Xd_pulse(:,i,j) = Xd_pulse(:,i,j) - smooth(Xd_pulse(:,i,j),length(Xd_pulse(:,i,j))/Nsm);
    end
end

for j = 1:a
    for i = 1:(logdata.data(1,iter_loc))
        NFFT = (2^8)*2^nextpow2(numel(time));
        nw = 1;
        XsegL = floor(length(time)/nw);
        welchwin = round(XsegL);
        NPSD = floor(NFFT/nw);%Prevents additional interpolation by pwelch
        noverlap = 0; %To generate an average spectrum of independent points use 0
        winfunc = rectwin(welchwin);
        %Find the window normalization factor for the peak amplitude
        freq = 0.005;
        Xsine = sin(2*pi*freq.*time);
        [Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
        winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
        [Xpsd(:,i,j),fpsd(:,i,j)] = pwelch(Xd_pulse(:,i,j),winfunc,noverlap,NPSD,Fs);
        fscale = 10^3;
        Xpsd(:,i,j) = Xpsd(:,i,j)./fscale;%Change units to (nm)^2/Hz
    end
end


% Find peaks
findpksyn=1;

fmin = 0.002;
fmax = 0.02;
elim = 0;               % eliminate 60, 120, and 180 Hz peaks?

if findpksyn ==1;
for j = 1:a
    for i = 1:(logdata.data(1,iter_loc))
        if elim == 1
            clear freqrange f60 f120 f180 freqrange2
            freqrange2 = find(fpsd(:,i,j) <= fmax & fpsd(:,i,j) >= fmin);
            f60 = find(fpsd(freqrange2,i,j) <= .061 & fpsd(freqrange2,i,j) >= .059); freqrange2(f60)=[];
            f120 = find(fpsd(freqrange2,i,j) <= .121 & fpsd(freqrange2,i,j) >= .119); freqrange2(f120)=[];
            f180 = find(fpsd(freqrange2,i,j) <= .181 & fpsd(freqrange2,i,j) >= .179); freqrange2(f180)=[];
            freqrange = (fpsd(:,i,j) <= fmax & fpsd(:,i,j) >= fmin & Xpsd(:,i,j) > circshift(Xpsd(:,i,j),[1 1]) & Xpsd(:,i,j) > circshift(Xpsd(:,i,j),[-1 1]));
            freqrange(:)=0;
            freqrange(freqrange2)=1;
        else
            clear freqrange
            freqrange= (fpsd(:,i,j) <= fmax & fpsd(:,i,j) >= fmin);
        end
Xpsdmaxima(:,i,j) = Xpsd(:,i,j).*(freqrange);

Xpsdpeak(i,j) = max(Xpsdmaxima(:,i,j));

XpsdmaximaMedian = median(Xpsdmaxima(find(Xpsdmaxima(:,i,j) >= 0.1*Xpsdpeak(i,j)),i,j));

%The main peak of bundle oscillations. Can be bigger than the local std.
XFTpeak(i,j) = (sqrt(fscale.*Xpsdpeak(i,j).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
fpsdpeak(i,j) = fpsd(findnearest(Xpsd(:,i,j), Xpsdpeak(i,j)));

    end
end
end
%}

%{
% Use multitaper method
nwmtm = 1.25;
for j = 1:a
    for i = 1:(logdata.data(1,iter_loc))
        [Xmtm(:,i,j) fmtm(:,i,j)] = pmtm(Xd_pulse(:,i,j),nwmtm,NFFT,Fs);
    end
end
%}

    
