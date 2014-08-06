function Tfind
%Experimental data. All points files%%%%%%%%%%%%%%%%%
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-04.01/Ear 1/Cell 4/20140804-cell4.mat');

%Operating points in ascending order
Fsort = sort(F_rand);
Fgrid = Fsort(diff(Fsort) ~= 0);
Fgrid(end+1) = max(F_rand);
ksort = sort(k_rand);
kgrid = ksort(diff(ksort) ~= 0);
kgrid(end+1) = max(k_rand);

Fgridrev = sort(Fgrid,'descend');

sizeXd = size(Xd);
Np = sizeXd(2);
Nt = sizeXd(3);

Tstartend(1:2,1:Np,1:Nt) = zeros(2,Np,Nt);

for Findex = 1:length(Fgrid)
for kindex = 1:length(kgrid)

Npt = find(k_rand == kgrid(kindex) & F_rand == Fgridrev(Findex));
if rem(Npt,Np) ~= 0;
    Npulse = rem(Npt,Np);
else
    Npulse = Np;
end  
Ntrial = (Npt - Npulse)/Np + 1;
Npt = (Ntrial-1)*Np+Npulse;

X = Xd(:,Npulse,Ntrial);

deltat = 0.1;%ms
tvec = 0:deltat:(length(X)-1)*deltat;

ssstartpt = 1;
ssendpt = length(X);

%Minimum window length allowed
Tmin  = min([3000 tvec(ssendpt)]); %ms

if tvec(ssendpt) - tvec(ssstartpt) >= Tmin
tmax = tvec(ssendpt);
tmin = tvec(ssstartpt);
else
tmax = tvec(ssendpt);
tmin = tmax - Tmin;
end

XpsdpeakMax = 0;
ratio = 0;
tstart = 0;
tend = 8400;
deltaT = 200;

tstartnew = tmin;
tendnew = tmax;

winsearch = 1;
if winsearch == 1
for m = 20:-1:-20
for l = -20:20
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

Tstartend(1,Npulse,Ntrial) = tstartnew;
Tstartend(2,Npulse,Ntrial) = tendnew;

end
end

%%%%%%%%%%%Save the time limits%%%%%%%%%%%
timefile = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-04.01/Ear 1/Cell 4/Tstartend.mat';
save(timefile, 'Tstartend');
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
fmin = 0.001;
%CHOICE
fmax = 0.1;

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
Xpsdmaxima = Xpsd.*(fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));

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
