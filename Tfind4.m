function Tfind4
%Experimental data. All points files%%%%%%%%%%%%%%%%%
%load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/20130908-cell15-2.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-05.01/Ear 1/Cell 11/20140805-cell11.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/Gentamicin/2014-08-05.01/Ear 1/Cell 4/20140805-cell4.mat')
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-12-18.01/Ear 1/Cell 1/Extracted Data.mat');



XdL = length(Xd);


    Tstartend(1:2,1:XdL) = zeros(2,XdL);
for Findex = 1:XdL



X = Xd{Findex};
X = X-smooth(X,length(X)/10);

deltat = 0.001;%ms
tvec = 0:deltat:(length(X)-1)*deltat;

ssstartpt = 1;
ssendpt = length(X);

%Minimum window length allowed
Tmin  = min([1500 tvec(ssendpt)]); %ms

if tvec(ssendpt) - tvec(ssstartpt) >= Tmin
tmax = tvec(ssendpt);
tmin = tvec(ssstartpt);
else
tmax = tvec(ssendpt);
tmin = tmax - Tmin;
end

XpsdpeakMax = 0;
ratio = 0;
tstart = 5000;
tend = 10000;
deltaT = 100;

tstartnew = tmin;
tendnew = tmax;

winsearch = 1;
if winsearch == 1
for m = 50:-1:-50
for l = -50:50
    if tstart+l*deltaT >= tmin && tend+m*deltaT <= tmax && Tmin <= tend+m*deltaT - (tstart+l*deltaT)
    [fpsdpeak, Xpsdpeak, XpsdmaximaMedian] = ssPSD(X,tstart+l*deltaT,tend+m*deltaT);
    if Xpsdpeak/XpsdmaximaMedian > ratio
        ratio = Xpsdpeak/XpsdmaximaMedian;
        XpsdpeakMax = Xpsdpeak;
        tstartnew = tstart + l*deltaT;
        tendnew = tend + m*deltaT;
        %'choice1'
    elseif Xpsdpeak/XpsdmaximaMedian == ratio && Xpsdpeak > XpsdpeakMax
        XpsdpeakMax = Xpsdpeak;
        tstartnew = tstart + l*deltaT;
        tendnew = tend + m*deltaT;
        %'choice2'
    elseif Xpsdpeak/XpsdmaximaMedian == ratio && Xpsdpeak == XpsdpeakMax && (tstartnew > tstart + l*deltaT || tendnew < tend + m*deltaT)
        tstartnew = tstart + l*deltaT;
        tendnew = tend + m*deltaT;
        %'choice3'
    end
    end
end
end
end

Tstartend(1,Findex) = tstartnew;
Tstartend(2,Findex) = tendnew;
[fpsdpeak(Findex), Xpsdpeak(Findex), XpsdmaximaMedian(Findex)] = ssPSD(X,tstartnew,tendnew);

if Tstartend(1,Findex) ~= Tstartend(2,Findex)
    ssstartpt = find(abs(tvec-Tstartend(1,Findex))==min(abs(tvec-Tstartend(1,Findex))));
    ssendpt = find(abs(tvec-Tstartend(2,Findex))==min(abs(tvec-Tstartend(2,Findex))));
    else
    ssstartpt = 1;
    ssendpt = round(length(Xd)/2);
    ssstartpt = 1;
    end
%Remove the mean
X = Xd{Findex}(ssstartpt:ssendpt)-mean(Xd{Findex}(ssstartpt:ssendpt));
tvec = tvec(ssstartpt:ssendpt);


%%%%%%%%%%%%%%%%PSD%%%%%%%%%%%%%%%%%
NFFT = (2^4)*2^nextpow2(numel(tvec));
%NFFT = numel(tvec);
nw = 1;
XsegL = floor(length(tvec)/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);%Prevents additional interpolation by pwelch
noverlap = 0; %To generate an average spectrum of independent points use 0
winfunc = hamming(welchwin);
%Find the window normalization factor for the peak amplitude
freq = 0.005;
Xsine = sin(2*pi*freq.*tvec);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
[Xpsd,fpsd] = pwelch(X,winfunc,noverlap,NPSD,Fs);

%Multitaper method distorts the spectrum greatly
%nw = 4;%Time-halfbandwidth product, 2*nw-1 tapers used, min(nw) = 1.25
%bwidth = 2*nw/tvec(end)
%[Xpsd,fpsd] = pmtm(X,nw,NPSD,Fs);

%Rescale the PSD
fscale = 10^3;
Xpsd = Xpsd./fscale;%Change units to (nm)^2/Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Same number of elements as Xpsd
% Eliminate 60,120,180 ± 1 Hz
elim=1;
fmax = 15;
fmin = 0.005;
if elim == 1
    clear freqrange f60 f120 f180 freqrange2
    freqrange2 = find(fpsd <= fmax & fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));
    f60 = find(fpsd(freqrange2) <= .061 & fpsd(freqrange2) >= .059); freqrange2(f60)=[];
    f120 = find(fpsd(freqrange2) <= .121 & fpsd(freqrange2) >= .119); freqrange2(f120)=[];
    f180 = find(fpsd(freqrange2) <= .181 & fpsd(freqrange2) >= .179); freqrange2(f180)=[];
    freqrange = (fpsd <= fmax & fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));
    freqrange(:)=0;
    freqrange(freqrange2)=1;
else
    clear freqrange
    freqrange= (fpsd <= fmax & fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));
end

Xpsdmaxima = Xpsd.*(freqrange);

Xpsdpeak = max(Xpsdmaxima);

XpsdmaximaMedian = median(Xpsdmaxima(find(Xpsdmaxima >= 0.1*Xpsdpeak)));

%The main peak of bundle oscillations. Can be bigger than the local std.
XFTpeak1(Findex) = (sqrt(fscale.*Xpsdpeak.*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
fpsdpeak1(Findex) = fpsd(find(Xpsd == Xpsdpeak));


end


%%%%%%%%%%%Save the time limits%%%%%%%%%%%
timefile = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-12-18.01/Ear 1/Cell 1/tstartendPSD.mat';
save(timefile, 'Tstartend','fpsdpeak','XFTpeak1','fpsdpeak1','tvec');
display('saving...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function [fpsdpeak, Xpsdpeak, XpsdmaximaMedian] = ssPSD(X,tmanstart,tmanend)
    
deltat = 0.1;%ms
Fs = 1/deltat;%kHz
tvec = 0:deltat:(length(X)-1)*deltat;

ssstartpt = find(abs(tvec-tmanstart)==min(abs(tvec-tmanstart)));
ssendpt = find(abs(tvec-tmanend)==min(abs(tvec-tmanend))); 

%Change X and tvec here
tvec = tvec(ssstartpt:ssendpt);
X = X(ssstartpt:ssendpt) - mean(X(ssstartpt:ssendpt));

%%%%%%%%%%%PSD Parameters%%%%%%%%%%%%%%%%
%fmin = Fs/length(Xd);%Actual frequency resolution if entire time trace is used
%CHOICE
%fmin = 0.002;
% CHOICE
fmin = 0.002;%kHz
%CHOICE
fmax = 100;

%%%%%%%%%%%%%%%%Windowed PSD%%%%%%%%%%%%%%%%%
NFFT = (2^4)*2^nextpow2(numel(tvec));
%NFFT = numel(tvec);

nw = 1;
XsegL = floor(length(tvec)/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);%Prevents additional interpolation by pwelch
noverlap = 0; %To generate an average spectrum of independent points use 0

%Changes spectrum as time window changes if frequency content is changing
%with time as the windows weight the middle of the trace more than the ends
%Should include parts of traces with less probable dynamics at ends

winfunc = hamming(welchwin,'periodic');
%4/ttot main lobe width, 8*10^-4 side lobe suppression
%0.5*ttot weighted above 50%
%=> 1 Hz res requires at least 4s

%Find the window normalization factor for the peak amplitude
freq = 0.005;
Xsine = sin(2*pi*freq.*tvec);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;

[Xpsd,fpsd] = pwelch(X,winfunc,noverlap,NPSD,Fs);

%Multitaper method distorts the spectrum greatly
%nw = 4;%Time-halfbandwidth product, 2*nw-1 tapers used, min(nw) = 1.25
%bwidth = 2*nw/tvec(end)
%[Xpsd,fpsd] = pmtm(X,nw,NPSD,Fs);

%Rescale the PSD
fscale = 10^3;
Xpsd = Xpsd./fscale;%Change units to (nm)^2/Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Same number of elements as Xpsd
% Eliminate 60,120,180 ± 1 Hz
elim=1;
if elim == 1
    clear freqrange f60 f120 f180 freqrange2
    freqrange2 = find(fpsd <= fmax & fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));
    f60 = find(fpsd(freqrange2) <= .061 & fpsd(freqrange2) >= .059); freqrange2(f60)=[];
    f120 = find(fpsd(freqrange2) <= .121 & fpsd(freqrange2) >= .119); freqrange2(f120)=[];
    f180 = find(fpsd(freqrange2) <= .181 & fpsd(freqrange2) >= .179); freqrange2(f180)=[];
    freqrange = (fpsd <= fmax & fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));
    freqrange(:)=0;
    freqrange(freqrange2)=1;
else
    clear freqrange
    freqrange= (fpsd <= fmax & fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));
end

%Same number of elements as Xpsd
Xpsdmaxima = Xpsd.*(freqrange);

Xpsdpeak = max(Xpsdmaxima);

%The main peak of bundle oscillations. Can be bigger than the local std.
XFTpeak = (sqrt(fscale.*Xpsdpeak.*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
fpsdpeak = fpsd(find(Xpsd == Xpsdpeak));

%Redefine Xpsdmaxima to include all the frequencies to choose a time segment
%with minimal drift
Xpsdmaxima = Xpsd.*(Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));

%The largest maxima excluding the peak
Xpsdmaxima(find(Xpsdmaxima==Xpsdpeak))=[];
%The median of the largest maxima
XpsdmaximaMedian = median(Xpsdmaxima(find(Xpsdmaxima >= 0.1*Xpsdpeak)));

end
