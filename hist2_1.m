clear all
close all

load('/Users/joshsalvi/Documents/Lab/Lab/Data Analysis/State Spaces/20130821/20130821-cell11.mat');

%time.data(:,1,1) is the time vec. Don't know what the other elements of
%time are but they're not time values.

AfterRamp = find(time.data(:,1,1) >= ramp);

%Use only second half of data to get steady state response
%CHOICE
AfterRamp(1) = round(length(time.data(:,1,1))/2);

BeforeRamp = find(time.data(:,1,1) <= time.data(end,1,1)-ramp);
timeR = time.data(AfterRamp(1):BeforeRamp(end),1,1);
clear time

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

modality = zeros(length(Fgrid),length(kgrid));

Fcpt = zeros(Np,Nt);
kpt = zeros(Np,Nt);

%Set pl = 1 to output plots
pl = 0;
for Npulse = 1:Np
    for Ntrial = 1:Nt
XdR = Xd(AfterRamp(1):BeforeRamp(end),Npulse,Ntrial);
if mean(XdR)==0
    break
end
%Modality of displacement histogram depends strongly upon how the data is
%detrended and upon whether the steady state has been reached.
%CHOICE
XdRS = XdR - smooth(XdR,round(length(XdR)/7),'sgolay',1);
XdRS = smooth(XdRS,Fs/9000,'sgolay',1);

%Bin width using the Freedman/Diaconis rule
bw = 2*iqr(XdRS)/length(XdRS)^(1/3);
Nbins = round((max(XdRS) - min(XdRS))/bw);

%if bw == 0
%    break;
%end
%Must redifine bin width such that Nbins*BinSize = max - min
BinSize = (max(XdRS) - min(XdRS))/Nbins;
[binFreq,binPos]=hist(XdRS,Nbins);
if pl == 1
    figure
    bar(binPos,binFreq)
end
%CHOICE
ws = round(Nbins/9); %Must be odd
if mod(ws,2) == 0
    ws = ws + 1;
end
%Changing the order of the polynomial changes the first derivative
%Looking for a cubic at most in the local distribution shape
[b,g] = sgolay(3,ws);   % Calculate S-G coefficients

dx = binPos(2) - binPos(1);
HalfWin  = ((ws+1)/2) -1;
SG0 = zeros(1,length((ws+1)/2:length(binFreq)-(ws+1)/2));
SG1 = zeros(1,length((ws+1)/2:length(binFreq)-(ws+1)/2));
SG2 = zeros(1,length((ws+1)/2:length(binFreq)-(ws+1)/2));
binPosS = zeros(1,length((ws+1)/2:length(binFreq)-(ws+1)/2));
for n = (ws+1)/2:length(binFreq)-(ws+1)/2,
  SG0(n-(ws+1)/2+1) = dot(g(:,1), binFreq(n - HalfWin: n + HalfWin));
  SG1(n-(ws+1)/2+1) = dot(g(:,2), binFreq(n - HalfWin: n + HalfWin));
  SG2(n-(ws+1)/2+1) = 2*dot(g(:,3)', binFreq(n - HalfWin: n + HalfWin))';
  binPosS(n-(ws+1)/2+1) = binPos(n);
end
if pl == 1
    hold on
    plot(binPosS,SG0,'g')
end

SG1 = SG1/dx;         % Turn differential into derivative
SG2 = SG2/dx^2;         % Turn differential into derivative

%Distribution must have positive slopes at beginning and end
binPosS(find(abs(diff(sign(SG1)))==2));
MaxNumber = (length(find(abs(diff(sign(SG1)))==2))+1)/2;
Npt = (Ntrial-1)*Np+Npulse;
modality_rand(Npt) = MaxNumber > 1;
modality(find(Fgrid == F_rand(Npt)),find(kgrid == k_rand(Npt))) = (MaxNumber > 1) + 1;
mod2_rand(Npt)=MaxNumber;
mod2(find(Fgrid == F_rand(Npt)),find(kgrid == k_rand(Npt))) = (MaxNumber > 2) + 1;

if pl == 1
    text(0.6*max(binPos),0.8*max(binFreq),[num2str(MaxNumber) ' max'],'HorizontalAlignment','left')
    text(0.6*max(binPos),0.75*max(binFreq),['Fc: ' num2str(Fcpt(Npulse,Ntrial)) ' pN'],'HorizontalAlignment','left')
    text(0.6*max(binPos),0.7*max(binFreq),['k: ' num2str(kpt(Npulse,Ntrial))  ' \muN\cdotm^{-1}'],'HorizontalAlignment','left')
end
    end
end

%Create grid so that data squares are centered correctly.
dFgrid = diff(Fgrid);
dFgrid(end+1) = dFgrid(end);
Fgrid = Fgrid-dFgrid/2;
Fgrid(end+1) = Fgrid(end) + dFgrid(end);
modality(end+1,:)=zeros(1,length(kgrid));
mod2(end+1,:)=zeros(1,length(kgrid));


dkgrid = diff(kgrid);
dkgrid(end+1) = dkgrid(end);
kgrid = kgrid-dkgrid/2;
kgrid(end+1) = kgrid(end) + dkgrid(end);
modality(:,end+1)=zeros(length(Fgrid),1);
mod2(:,end+1)=zeros(length(Fgrid),1);

kgrid=min(kgrid):(max(kgrid)-min(kgrid))/(length(kgrid)-1):max(kgrid);

figure
pcolor(kgrid*10^6,Fgrid*10^12,modality)
axis square;
%colormap(gray(2))
title('Modality'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');
colorbar

figure
pcolor(kgrid*10^6,Fgrid*10^12,mod2)
axis square;
%colormap(gray(2))
title('Multimodal'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');
colorbar

%figure
%scatter3(k_rand*10^6,F_rand*10^12,modality_rand,100*ones(size(modality_rand)),modality_rand,'filled')
%view([0 90])
%title('Modality'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');
