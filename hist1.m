clear all
%close all

load('/Users/joshsalvi/Documents/Lab/Lab/Data Analysis/State Spaces/20130821/20130821-cell11.mat');

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

kpt = kgrid(5);
Fpt = Fgrid (5);
kpt*10^6
Fpt*10^12
Nptlist = find(k_rand == kpt & F_rand == Fpt);
Npt = Nptlist(1);
Npulse = mod(Npt-1,Np)+1
Ntrial = ceil(Npt/Np)

%plot control parameters and responses on a normalized scale
plot(time.data(:,1,1),Xc(:,Npulse,Ntrial)/max(abs(Xc(:,Npulse,Ntrial))),'k')
hold on
%plot(time.data(:,1,1),Gain(:,Npulse,Ntrial)/max(abs(Gain(:,Npulse,Ntrial))),'g')
plot(time.data(:,1,1),Xd(:,Npulse,Ntrial)/max(abs(Xd(:,Npulse,Ntrial))),'r')
%plot(time.data(:,1,1),Delta(:,Npulse,Ntrial)/max(abs(Delta(:,Npulse,Ntrial))),'b')

AfterRamp = find(time.data(:,1,1) >= ramp*100);

%Use only end of data to get steady state response
AfterRamp(1) = round(0.5*length(time.data(:,1,1)));

BeforeRamp = find(time.data(:,1,1) <= time.data(end,1,1)-ramp*100);

XdR = Xd(AfterRamp(1):BeforeRamp(end),Npulse,Ntrial);
timeR = time.data(AfterRamp(1):BeforeRamp(end),1,1);

%Modality of displacement histogram depends strongly upon how the data is
%detrended an upon whether the steady state has been reached.
XdRS = XdR - smooth(XdR,round(length(XdR)/5),'sgolay',1);
%The mean should be near zero
mean(XdRS);
figure
plot(timeR,XdRS,'r')
XdRSS = smooth(XdRS,Fs/2000,'sgolay',1);
hold on
plot(timeR,XdRSS,'k')
XdRS = XdRSS;

%Bin width using the Freedman?Diaconis rule
bw = 2*iqr(XdRS)/length(XdRS)^(1/3);
Nbins = round((max(XdRS) - min(XdRS))/bw);
%Must redifine bin width such that Nbins*BinSize = max - min
BinSize = (max(XdRS) - min(XdRS))/Nbins;
[binFreq,binPos]=hist(XdRS,Nbins);
figure
bar(binPos,binFreq)
ws = round(Nbins/6); %Must be odd
if mod(ws,2) == 0
    ws = ws + 1;
end
%binFreqS = smooth(binFreq,ws,'sgolay',1);
%hold on
%plot(binPos,binFreqS,'r')

%Changing the order of the polynomial changes the first derivative
%Looking for a cubic at most in the local distribution shape
[b,g] = sgolay(3,ws);   % Calculate S-G coefficients

dx = binPos(2) - binPos(1);
HalfWin  = ((ws+1)/2) -1;
SG0 = zeros(1,length((ws+1)/2:length(binFreq)-(ws+1)/2));
SG1 = zeros(1,length((ws+1)/2:length(binFreq)-(ws+1)/2));
SG2 = zeros(1,length((ws+1)/2:length(binFreq)-(ws+1)/2));
for n = (ws+1)/2:length(binFreq)-(ws+1)/2,
  SG0(n-(ws+1)/2+1) = dot(g(:,1), binFreq(n - HalfWin: n + HalfWin));
  SG1(n-(ws+1)/2+1) = dot(g(:,2), binFreq(n - HalfWin: n + HalfWin));
  SG2(n-(ws+1)/2+1) = 2*dot(g(:,3)', binFreq(n - HalfWin: n + HalfWin))';
  binPosS(n-(ws+1)/2+1) = binPos(n);
end
hold on
plot(binPosS,SG0,'g')

SG1 = SG1/dx;         % Turn differential into derivative
SG2 = SG2/dx^2;         % Turn differential into derivative
figure
plot(binPosS,SG0,'g')
hold on
plot(binPosS,SG1,'r')
plot(binPosS,SG2,'b')
line([binPosS(1) binPosS(end)],[0 0],'LineStyle','-')

%Distribution must have positive slopes at beginning and end
binPosS(find(abs(diff(sign(SG1)))==2))
MaxNumber = (length(find(abs(diff(sign(SG1)))==2))+1)/2
text(0.8*max(binPos),0.8*max(binFreq),[num2str(MaxNumber) ' max'],'HorizontalAlignment','left')

h = 1;
if h == 1
XdRSA = hilbert(XdRS);
minR = min(real(XdRSA));
maxR = max(real(XdRSA));
minI = min(imag(XdRSA));
maxI = max(imag(XdRSA));
NRbins = round(maxR - minR)*1;
NIbins = round(maxI - minI)*1;
BinSizeR = ((maxR - minR))/NRbins;
BinSizeI = ((maxI - minI))/NIbins;
AHist = zeros(NRbins,NIbins);
Rpos = (minR+BinSizeR/2):BinSizeR:(maxR-BinSizeR/2);
Ipos = (minI+BinSizeI/2):BinSizeI:(maxI-BinSizeI/2);

for R = 1:NRbins
    for I = 1:NIbins
    AHist(R,I) = sum(real(XdRSA)<(Rpos(R)+BinSizeR/2) & real(XdRSA)>=(Rpos(R)-BinSizeR/2) & imag(XdRSA)<(Ipos(I)+BinSizeI/2) & imag(XdRSA)>=(Ipos(I)-BinSizeI/2));
    end
end
figure
surf(Ipos,Rpos,AHist)

TestT = 100; %Period in ms
w = 2*pi/TestT; %Time in ms
TestTime = 0:TestT/1000:10*TestT;

%COSINE
%TestSigR = cos(w.*TestTime) + 0.1.*randn(length(TestTime),1)';

%SQUARE WAVE
duty = 50; %Percentange of period signal is positive
%TestSigR = square(w.*TestTime,duty) + 0.1.*randn(length(TestTime),1)';

%RELAXATION OSCILLATION
limcycnoise;
TestSigR = fliplr(y_amp(:,1,1)');
TestTime = fliplr(t_tot(:,:,1)');

TestSig = hilbert(TestSigR);
minR = min(real(TestSig));
maxR = max(real(TestSig));
minI = min(imag(TestSig));
maxI = max(imag(TestSig));
NRbins = round(maxR - minR)*10;
NIbins = round(maxI - minI)*10;
BinSizeR = ((maxR - minR))/NRbins;
BinSizeI = ((maxI - minI))/NIbins;
AHist = zeros(NRbins,NIbins);
Rpos = (minR+BinSizeR/2):BinSizeR:(maxR-BinSizeR/2);
Ipos = (minI+BinSizeI/2):BinSizeI:(maxI-BinSizeI/2);

for R = 1:NRbins
    for I = 1:NIbins
    AHist(R,I) = sum(real(TestSig)<(Rpos(R)+BinSizeR/2) & real(TestSig)>=(Rpos(R)-BinSizeR/2) & imag(TestSig)<(Ipos(I)+BinSizeI/2) & imag(TestSig)>=(Ipos(I)-BinSizeI/2));
    end
end
figure
surf(Ipos,Rpos,AHist)
figure
plot(TestTime,real(TestSig))
end
