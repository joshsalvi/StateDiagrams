clear all
%close all

load('/Users/joshsalvi/Documents/Lab/Lab/Data Analysis/State Spaces/20130507/20130507-cell3.mat');

%time.data(:,1,1) is the time vec. Don't know what the other elements of
%time are but they're not time values.
sizeXd = size(Xd);
Np = sizeXd(2);
Nt = sizeXd(3);

%Three points are duplicates at (0,0) 
Fsort = sort(F_rand);
Fgrid = Fsort(diff(Fsort) ~= 0);
Fgrid(end+1) = max(F_rand);
ksort = sort(k_rand);
kgrid = ksort(diff(ksort) ~= 0);
kgrid(end+1) = max(k_rand);

kpt = kgrid(4);
Fpt = Fgrid (4);
kpt*10^6
Fpt*10^12
Nptlist = find(k_rand == kpt & F_rand == Fpt);
Npt = Nptlist(1);
Npulse = mod(Npt-1,Np)+1
Ntrial = ceil(Npt/Np)

%plot control parameters and responses on a normalized scale
figure
plot(time.data(:,1,1),Xc(:,Npulse,Ntrial)/max(abs(Xc(:,Npulse,Ntrial))),'k')
hold on
%plot(time.data(:,1,1),Gain(:,Npulse,Ntrial)/max(abs(Gain(:,Npulse,Ntrial))),'g')
plot(time.data(:,1,1),Xd(:,Npulse,Ntrial)/max(abs(Xd(:,Npulse,Ntrial))),'r')
%plot(time.data(:,1,1),Delta(:,Npulse,Ntrial)/max(abs(Delta(:,Npulse,Ntrial))),'b')

AfterRamp = find(time.data(:,1,1) >= ramp);
BeforeRamp = find(time.data(:,1,1) <= time.data(end,1,1)-ramp);
AfterRamp(1) = round(0.5*length(time.data(:,1,1)));

XdR = Xd(AfterRamp(1):BeforeRamp(end),Npulse,Ntrial);
timeR = time.data(AfterRamp(1):BeforeRamp(end),1,1);

test = 0;
if test == 1
%Test functions
%Sinusoidal wave
TestT = 100; %Period in ms
w = 2*pi/TestT; %Time in ms
%XdR = cos(w.*timeR) + 0.1.*randn(length(timeR),1);
%XdR = 0.1.*randn(length(timeR),1);
%Square wave. %Harmonic amplitude decreases as 1/f for 50% square wave but more slowly for asymm wave!
duty = 10; %Percentange of period signal is positive
XdR = square(w.*timeR,duty) + 0.1.*randn(length(timeR),1);
%Rectangular pulse. FT is sinc(f) with zeros at min(1/(pulse duration), 1/(total time - pulse duration)).
%pulse_start = 5000;
%pulse_end = length(timeR) - 5000;
%XdR = zeros(length(timeR),1);
%XdR(pulse_start:pulse_end) = ones(length(pulse_start:pulse_end),1);
end

%detrended data
XdRS = XdR - smooth(XdR,round(length(XdR)/10),'sgolay',1);

%Smooth the time series without detrending
XdRSS = smooth(XdR,round(0.01*length(XdR)),'sgolay',1);

figure
%plot(timeR,XdRS,'r')
hold on
plot(timeR,XdR,'g')
%plot(timeR,XdRSS,'k')
%XdR = XdRSS;

%Fourier transform
NFFT = (2^4)*2^nextpow2(numel(timeR));
%NFFT = numel(timeR);
f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT
actual_delta_f = Fs/numel(timeR)
%double sided
XdRFT = fft(XdR,NFFT)/numel(timeR);
XdRSFT = fft(XdRS,NFFT)/numel(timeR);
XdRSSFT = fft(XdRSS,NFFT)/numel(timeR);
%single sided
XdRFTss = 2*XdRFT(1:length(f));
XdRFTss(1) = XdRFTss(1)/sqrt(2); %0 freq. component increases by sqrt(2)
XdRFTss(end) = XdRFTss(end)/sqrt(2); %last freq. component increases by sqrt(2)

XdRSFTss = 2*XdRSFT(1:length(f));
XdRSFTss(1) = XdRSFTss(1)/sqrt(2); %0 freq. freq. component increases by sqrt(2)
XdRSFTss(end) = XdRSFTss(end)/sqrt(2); %last freq. component increases by sqrt(2)

XdRSSFTss = 2*XdRSSFT(1:length(f));
XdRSSFTss(1) = XdRSSFTss(1)/sqrt(2); %0 freq. freq. component increases by sqrt(2)
XdRSSFTss(end) = XdRSSFTss(end)/sqrt(2); %last freq. component increases by sqrt(2)

%Check to see what repeating the time series does to the spectrum
%Repeating the time series adds zeros between each pair of FT points.
%Don't add any new frequency components by repeating time series.
repeat = 0;
if repeat == 1
XdRrep = [XdR' XdR' XdR']';
timeRrep = [timeR' (timeR(end) + timeR - timeR(1))']';
timeRrep = [timeRrep' (timeRrep(end) + timeR - timeRrep(1))']';
figure
plot(timeRrep,XdRrep,'k')
hold on
plot(timeR,XdR,'r')
NFFTrep = numel(timeRrep);
frep = (Fs/2)*linspace(0,1,NFFTrep/2+1);
%double sided
XdRrepFT = fft(XdRrep,NFFTrep)/numel(timeRrep);
%single sided
XdRrepFTss = 2*XdRrepFT(1:length(frep));
XdRrepFTss(1) = XdRrepFTss(1)/sqrt(2); %0 freq. component increases by sqrt(2)
XdRrepFTss(end) = XdRrepFTss(end)/sqrt(2); %last freq. component increases by sqrt(2)
figure
plot(f,abs(XdRFTss),'red')
hold on
plot(frep,abs(XdRrepFTss),'k')
set(gca,'XScale','log')
set(gca,'YScale','log')
figure
plot(frep,abs(XdRrepFTss),'k')
set(gca,'XScale','log')
set(gca,'YScale','log')
return
end

%Detrending the time series changes the FT too much
figure
plot(f,abs(XdRFTss),'red')
hold on
%plot(f,abs(XdRSSFTss),'black')
plot(f,1./(f),'g')
set(gca,'XScale','log')
set(gca,'YScale','log')

%Find the peak in the spectrum for frequencies bigger than resolution
minf = 10*actual_delta_f;
XdRmodspec = abs(XdRFTss).*(f>minf)';
'FT peak'
f(find(XdRmodspec == max(XdRmodspec)))
abs(XdRFTss(find(XdRmodspec == max(XdRmodspec))))

fitted = 0;
if fitted == 1
%Fit the PSD to find the measurement filter%%%%%%%%%%%
s = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[3000 0]);
%s = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[3000]);

%Lorenzian PSD
fitF = fittype('a/(f0sq+x^2)','options',s);
%fitF = fittype('a/(10^5+x^2)','options',s);

%Fit the data using the fit options
[PSDfit,gof,output] = fit(abs(XdRFTss).^2,f',fitF)

figure
plot(f,abs(XdRFTss).^2,'red')
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
plot(PSDfit,'g')
plot(f,1./(f.^2),'k')
kB = 1.3806503*10^-23;
T = 295;
%f0 = k/(2*pi*gamma), a = kB*T/(pi^2*gamma)
PSDunits = 10^-18; %m^2/Hz
gamma = kB*T/(pi^2*PSDfit.a*PSDunits)
k = sqrt(PSDfit.f0sq)*2*pi*gamma
k_rand((Ntrial-1)*Np+Npulse)
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Smooth the FT
XdRFTabsS = smooth(abs(XdRFTss),10*round(actual_delta_f/effective_delta_f),'sgolay',1);
figure
plot(f,abs(XdRFTss),'red')
hold on
plot(f,XdRFTabsS,'black')
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Smoothed FT'); xlabel('Frequency (Hz)'); ylabel('X (nm)');
return

%Detrend the FT by subtracting the smoothed FT
XdRFTS = smooth(XdRFTss,200*round(actual_delta_f/effective_delta_f),'sgolay',1);
XdRFTDT = XdRFTss-XdRFTS;
figure
plot(f,abs(XdRFTss),'red')
hold on
plot(f,abs(XdRFTS),'black')
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Smoothed FT'); xlabel('Frequency (Hz)'); ylabel('X (nm)');
figure
plot(f,abs(XdRFTDT),'black')
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Detrended FT'); xlabel('Frequency (Hz)'); ylabel('X (nm)');

%Find the peak of the unmodified spectrum near the peak
%of the modified one
maxf = 10^2;
XdRmodspec = abs(XdRFTDT).*(f>minf)';
fpt = find(XdRmodspec == max(XdRmodspec));
fpoint_min = find( abs(f-(f(fpt)-actual_delta_f)) == min(abs(f-(f(fpt)-actual_delta_f))) );
fpoint_max = find( abs(f-(f(fpt)+actual_delta_f)) == min(abs(f-(f(fpt)+actual_delta_f))) );
'FT peak near FT detrended by f peak'
f(find(abs(XdRFTss) == max(abs(XdRFTss(fpoint_min:fpoint_max)))))
max(abs(XdRFTss(fpoint_min:fpoint_max)))
XdRFTmodspecSTD = std(abs(XdRmodspec(fpoint_min:fpoint_max)))
XdRFTmodspecM = mean(abs(XdRmodspec(fpoint_min:fpoint_max)))
max(abs(XdRmodspec(fpoint_min:fpoint_max)))-XdRFTmodspecM > 2*XdRFTmodspecSTD

%Find the peak of the unmodified spectrum near the peak
%of the modified one
XdRmodspec = abs(XdRSSFTss).*(f>minf)';
fpt = find(XdRmodspec == max(XdRmodspec));
fpoint_min = find( abs(f-(f(fpt)-actual_delta_f)) == min(abs(f-(f(fpt)-actual_delta_f))) );
fpoint_max = find( abs(f-(f(fpt)+actual_delta_f)) == min(abs(f-(f(fpt)+actual_delta_f))) );
'FT peak near FT*f peak of smoothed data'
f(find(abs(XdRFTss) == max(abs(XdRFTss(fpoint_min:fpoint_max)))))
max(abs(XdRFTss(fpoint_min:fpoint_max)))
XdRFTmodspecSTD = std(abs(XdRmodspec(fpoint_min:fpoint_max)))
XdRFTmodspecM = mean(abs(XdRmodspec(fpoint_min:fpoint_max)))
max(abs(XdRmodspec(fpoint_min:fpoint_max)))-XdRFTmodspecM > 2*XdRFTmodspecSTD

%Compare PSD to FT
%sqrt(PSD) same as single-sided FT 
%h = spectrum.periodogram('Rectangular');
%Hpsd=psd(h,XdR,'NFFT',NFFT,'Fs',Fs,'ConfLevel',0.95);
%figure
%plot(Hpsd.Frequencies,sqrt(Hpsd.Data),'red');
%hold on
%plot(f,abs(XdRFTss),'black')
%set(gca,'XScale','log')
%set(gca,'YScale','log')
%figure
%plot(f,abs(XdRFTss)./sqrt(Hpsd.Data),'black')

%Rectangular window has best resolution but largest scalloping loss due
%to sampling away from spectrum peak. Scalloping can be reduced by zero
%padding or using lower resolution windows. Rectangular window is best
%option for detecting sinusoidal signal from white noise.
%Use Welch's method to find the PSD
window=round(length(XdR)/4); %Important to get window size right
poverlap=0;%percentange overlap. 0% is Barlett's method
%Should only window before zero padding, not after!
h = spectrum.welch('Rectangular',window,poverlap);
Hpsd=psd(h,XdR,'NFFT',NFFT,'Fs',Fs,'ConfLevel',0.95);
figure
plot(Hpsd.Frequencies,sqrt(Hpsd.Data),'red');
hold on
plot(Hpsd.Frequencies,1./Hpsd.Frequencies,'g');
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Welch FT'); xlabel('Frequency (Hz)'); ylabel('X (nm)');

XdRmodspec = sqrt(Hpsd.Data).*(f>minf)';
'FT peak using Welch'
f(find(XdRmodspec == max(XdRmodspec)))
sqrt(Hpsd.Data(find(XdRmodspec == max(XdRmodspec))))

%RMS amplitude of detrended data
XdRSm = mean(XdRS);
XdRSrms = sqrt(sum((XdRS-XdRSm).^2)/length(XdRS))
