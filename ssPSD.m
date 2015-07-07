function [fpsdpeak, Xpsdpeak, XpsdmaximaMedian] = ssPSD(X,tmanstart,tmanend,Fs)
    
deltat = 1/(Fs*1e-3);%ms
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
