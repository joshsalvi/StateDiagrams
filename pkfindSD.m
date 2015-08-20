function pkfindSD()

%Experimental data. All points files%%%%%%%%%%%%%%%%%
%clear;close all;clc;
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Extracted Data-good2.mat')
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 6/Tstartend4-good2.mat');
figure;
fishfigs = 1;
if fishfigs == 1
%load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/Gentamicin/2014-08-05.01/Ear 1/Cell 8/Modality-foranalysis.mat');
modfile = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Modality2sec2Hzmin-good-dipstat-filtered.mat';
load(modfile);
load('/Users/joshsalvi/GitHub/StateDiagrams/customcolormaps-redblue.mat');
end
%%%%%%%%%%%%%%%% CHOOSE THRESHOLD %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 threshold = 25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%}
%{
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/Gentamicin/2014-08-05.01/Ear 1/Cell 7/20130805-cell7.mat');
fishfigs = 1;
if fishfigs == 1
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/4-Modality2sec1Hzmin.mat');
load('/Users/joshsalvi/GitHub/StateDiagrams/customcolormaps-redblue.mat');
end
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/4-Tstartend2sec1Hzmin.mat');
%}
%Operating points in ascending order
Fsort = sort(F_rand);
Fgrid = Fsort(diff(Fsort) ~= 0);
Fgrid(end+1) = max(F_rand);
ksort = sort(k_rand);
kgrid = ksort(diff(ksort) ~= 0);
kgrid(end+1) = max(k_rand);

sizeXd = size(Xd);
Np = sizeXd(2);
Nt = sizeXd(3);

deltat = 1/(Fs*1e-3); %ms
tvec = 0:deltat:(length(Xd)-1)*deltat;
tmin = 0;%ms
tmax = tvec(end);

%%%%%%%%%%%%%%%%%%%%Steady-state%%%%%%%%%%%
%Check that the local mean and std have converged

%Remove the drift from each time trace
fdrift = 0.001;%Hz
%ws must be odd
ws = (1/fdrift)/deltat;
if rem(floor(ws),2) == 1 %odd
    ws = floor(ws);
else %even
    ws = floor(ws)+1; 
end
%Find the running mean and std for a time greater than the longest period (1/fdrift)
tss = 1.5*(1/fdrift);
ssws = round(tss/deltat);
nstdss = 2;
stdfrac = 0.1;%Get to 90% of std ss
%Discount end discontinuities
ndev = 2;

%%%%%%%%%%%PSD Parameters%%%%%%%%%%%%%%%%
%fmin = Fs/length(Xd);%Actual frequency resolution if entire time trace is used
fmin = 2;
fmax = 50;
elim = 0;               % eliminate 60, 120, and 180 Hz peaks?
Xpsdminlim = 10^-1;
Xpsdmaxlim = 10^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
set(gca, 'LooseInset', get(gca, 'TightInset'));
Pwidth = Ssize(3);
Pheight = Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

Fgridrev = sort(Fgrid,'descend');

%Frequency grid 
freqgrid = zeros(length(Fgrid),length(kgrid));
%Amplitude grid 
freqpk = zeros(length(Fgrid),length(kgrid));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Grid Loops%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for Findex = 1:length(Fgrid)
%for kindex = 1:length(kgrid)
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

tvec = 0:deltat:(length(Xd)-1)*deltat;

    if Tstartend(1,Npulse,Ntrial) ~= Tstartend(2,Npulse,Ntrial)
    ssstartpt = find(abs(tvec-Tstartend(1,Npulse,Ntrial))==min(abs(tvec-Tstartend(1,Npulse,Ntrial))));
    ssendpt = find(abs(tvec-Tstartend(2,Npulse,Ntrial))==min(abs(tvec-Tstartend(2,Npulse,Ntrial))));
    else
    ssstartpt = 1;
    ssendpt = round(length(Xd)/2);
    ssstartpt = 1;
    end

 %%%%%%%%Steady State%%%%%%%%%%%%
    Xnodrift = Xd(1:ssendpt,Npulse,Ntrial) - smooth(Xd(1:ssendpt,Npulse,Ntrial),ws,'sgolay',1);
    tvec = 0:deltat:(length(Xd)-1)*deltat;
if Mod(Npulse,Ntrial) == 1    
    % Find peaks
[pk{Findex,kindex} tr{Findex,kindex}] = PTDetect(Xnodrift,threshold);
freqpk(Findex,kindex) = length(pk{Findex,kindex})/(tvec(length(Xnodrift))/1e3);
freqtr(Findex,kindex) = length(tr{Findex,kindex})/(tvec(length(Xnodrift))/1e3);

for j = 2:length(pk{Findex,kindex})
    pkt(j-1) = tvec(pk{Findex,kindex}(j)) - tvec(pk{Findex,kindex}(j-1));
end
for j = 2:length(tr{Findex,kindex})
    trt(j-1) = tvec(tr{Findex,kindex}(j)) - tvec(tr{Findex,kindex}(j-1));
end

CVpk(Findex,kindex) = std(pkt/1e3)/mean(pkt/1e3)*2;
CVtr(Findex,kindex) = std(trt/1e3)/mean(trt/1e3)*2;

if length(tr{Findex,kindex}) >= length(pk{Findex,kindex})
for j = 1:length(pk{Findex,kindex})
    amplpktr{Findex,kindex}(j) = abs(Xnodrift(pk{Findex,kindex}(j)) - Xnodrift(tr{Findex,kindex}(j)));
end
else
    for j = 1:length(tr{Findex,kindex})
    amplpktr{Findex,kindex}(j) = abs(Xnodrift(pk{Findex,kindex}(j)) - Xnodrift(tr{Findex,kindex}(j)));
end
end
pktrampl(Findex,kindex) = mean(amplpktr{Findex,kindex});
else
    pktrampl(Findex,kindex) = 0;
    amplpktr{Findex,kindex} = NaN;
    freqpk(Findex,kindex) = 0;
    freqtr(Findex,kindex) = 0;
    [pk{Findex,kindex} tr{Findex,kindex}] = PTDetect(Xnodrift,threshold);
for j = 2:length(pk{Findex,kindex})
    pkt(j-1) = tvec(pk{Findex,kindex}(j)) - tvec(pk{Findex,kindex}(j-1));
end
for j = 2:length(tr{Findex,kindex})
    trt(j-1) = tvec(tr{Findex,kindex}(j)) - tvec(tr{Findex,kindex}(j-1));
end

CVpk(Findex,kindex) = std(pkt/1e3)/mean(pkt/1e3)*2;
CVtr(Findex,kindex) = std(trt/1e3)/mean(trt/1e3)*2;
end





end
end

if fishfigs == 0
return
end

%%%%%%%%%%%%%%%%%%%%Fish Figures%%%%%%%%%%%%%%

%Flip grid along the force direction
freqgrid = flipdim(freqgrid,1);
freqpk = flipdim(freqpk,1);
freqtr = flipdim(freqtr,1);
pktrampl = flipdim(pktrampl,1);
CVpk = flipdim(CVpk,1);
CVtr = flipdim(CVtr,1);

pktrampl = pktrampl./2;

%%%%%Frequency and amplitude correlations%%%%%%%%%
freqvec = reshape(freqpk,[length(freqpk(:,1))*length(freqpk(1,:)) 1]);
%Find the points that are within the fish
fishvecpts = find(freqvec ~=0);
freqvec = freqvec(fishvecpts);
amplvec = reshape(pktrampl,[length(pktrampl(:,1))*length(pktrampl(1,:)) 1]);
%Find the points that are within the fish
fishvecpts = find(amplvec ~=0);
amplvec = amplvec(fishvecpts);


Fgrid2 = Fgrid'*ones(1,length(kgrid));
kgrid2 = ones(length(kgrid),1)*Fgrid;
Fvec = reshape(Fgrid2,[length(Fgrid2(:,1))*length(kgrid2(:,1)) 1]);
Fvec = Fvec(fishvecpts);
kvec = reshape(kgrid2,[length(kgrid2(:,1))*length(Fgrid2(:,1)) 1]);
kvec = kvec(fishvecpts);
[rhofreqampl,prhofreqampl]=corr(freqvec,amplvec,'type','Spearman');
[rhofreqk,prhofreqk]=corr(freqvec,kvec,'type','Spearman');
[rhofreqF,prhofreqF]=corr(freqvec,Fvec,'type','Spearman');
[rhoamplk,prhoamplk]=corr(amplvec,kvec,'type','Spearman');
[rhoamplF,prhoamplF]=corr(amplvec,Fvec,'type','Spearman');

save(['/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/FreqAmplcorrelations2sec2Hzmin-good2-tfind2-pkdetect-threshold' num2str(threshold) '.mat'],'rhofreqampl','prhofreqampl','rhofreqk','prhofreqk','rhofreqF','prhofreqF','rhoamplk','prhoamplk','rhoamplF','prhoamplF','freqpk','freqtr','CVtr','CVpk','pktrampl','pk','tr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create grid so that data squares are centered correctly.
dFgrid = diff(Fgrid);
dFgrid(end+1) = dFgrid(end);
Fgridold = Fgrid;
Fgrid = Fgrid-dFgrid/2;
Fgrid(end+1) = Fgrid(end) + dFgrid(end);
freqgrid(end+1,:)=zeros(1,length(kgrid));
freqpk(end+1,:)=zeros(1,length(kgrid));
freqtr(end+1,:)=zeros(1,length(kgrid));
pktrampl(end+1,:)=zeros(1,length(kgrid));
CVpk(end+1,:)=zeros(1,length(kgrid));
CVtr(end+1,:)=zeros(1,length(kgrid));

dkgrid = diff(kgrid);
dkgrid(end+1) = dkgrid(end);
kgridold = kgrid;
kgrid = kgrid-dkgrid/2;
kgrid(end+1) = kgrid(end) + dkgrid(end);
freqgrid(:,end+1)=zeros(length(Fgrid),1);
freqpk(:,end+1)=zeros(length(Fgrid),1);
freqtr(:,end+1)=zeros(length(Fgrid),1);
pktrampl(:,end+1)=zeros(length(Fgrid),1);
CVtr(:,end+1)=zeros(length(Fgrid),1);
CVpk(:,end+1)=zeros(length(Fgrid),1);

kgrid=min(kgrid):(max(kgrid)-min(kgrid))/(length(kgrid)-1):max(kgrid);

FS = 24; %Font Size
FN = 'Arial'; %Font Name
LW = 2; %LineWidth
TL = 0.015; %TickLength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Frequency Figure%%%%%%%%%%%%%%
fh = figure;

set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
Pwidth = 0.45*Ssize(3);
Pheight = 1*Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

fh = pcolor(kgrid,Fgrid,freqpk);
colormap(blue1);

cd=get(fh,'cdata');
cd(find(freqpk==0))=NaN;
set(fh,'cdata',cd)

%The color range is rescaled to match the nonzero data range
caxis([min(min(freqpk(find(freqpk~=0)))) max(max(freqpk))])
axis square;
xlabel({'Stiffness (mN\cdotm^{-1})' ''},'FontSize',FS,'FontName',FN);
ylabel('Force (pN)','FontSize',FS,'FontName',FN);
kticks = ((kgrid+circshift(kgrid,[1 -1]))/2)*10^1;
kticks(end) = [];
Fticks = ((Fgrid+circshift(Fgrid,[1 -1]))/2)*10^1;
Fticks(end) = [];
%Set the labels at the correct positions and round the label values
set(gca,'xtick',kticks,'ytick',Fticks,'FontSize',FS,'FontName',FN)
set(gca,'xticklabel',round(kticks),'ytick',round(Fticks),'FontSize',FS,'FontName',FN)

set(gca,'TickDir','Out')
set(gca, 'TickLength', [TL TL]);
set(gca,'box','off')
set(gca,'LineWidth',LW)
set(fh,'LineWidth',LW)

cbh = colorbar('location','Southoutside');
xlabel(cbh,'Amplitude (nm)','FontSize',FS,'FontName',FN);
set(gca,'FontSize',FS,'FontName',FN)
set(cbh, 'TickLength', [TL TL]);
set(cbh,'TickDir','Out')
set(cbh,'box','off')
set(cbh,'LineWidth',LW)
set(gca,'LineWidth',LW)
set(cbh,'xtick',ceil(min(min(freqpk(find(freqpk~=0))))):ceil((floor(max(max(freqpk)))-ceil(min(min(freqpk(find(freqpk~=0))))))/10):floor(max(max(freqpk))))
set(cbh,'xticklabel',ceil(min(min(freqpk(find(freqpk~=0))))):ceil((floor(max(max(freqpk)))-ceil(min(min(freqpk(find(freqpk~=0))))))/10):floor(max(max(freqpk))))

set(gca, 'LooseInset', [0,0.2,0,0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Amplitude Figure%%%%%%%%%%%%%%
fh = figure;

set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
Pwidth = 0.45*Ssize(3);
Pheight = 1*Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

fh = pcolor(kgrid,Fgrid,pktrampl);
colormap(red1);

cd=get(fh,'cdata');
cd(find(pktrampl==0))=NaN;
set(fh,'cdata',cd)

%The color range is rescaled to match the nonzero data range
caxis([min(min(pktrampl(find(pktrampl~=0)))) max(max(pktrampl))])
axis square;
xlabel({'Stiffness (mN\cdotm^{-1})' ''},'FontSize',FS,'FontName',FN);
ylabel('Force (pN)','FontSize',FS,'FontName',FN);
kticks = ((kgrid+circshift(kgrid,[1 -1]))/2)*10^1;
kticks(end) = [];
Fticks = ((Fgrid+circshift(Fgrid,[1 -1]))/2)*10^1;
Fticks(end) = [];
%Set the labels at the correct positions and round the label values
set(gca,'xtick',kticks,'ytick',Fticks,'FontSize',FS,'FontName',FN)
set(gca,'xticklabel',round(kticks),'ytick',round(Fticks),'FontSize',FS,'FontName',FN)

set(gca,'TickDir','Out')
set(gca, 'TickLength', [TL TL]);
set(gca,'box','off')
set(gca,'LineWidth',LW)
set(fh,'LineWidth',LW)

cbh = colorbar('location','Southoutside');
xlabel(cbh,'Amplitude (nm)','FontSize',FS,'FontName',FN);
set(gca,'FontSize',FS,'FontName',FN)
set(cbh, 'TickLength', [TL TL]);
set(cbh,'TickDir','Out')
set(cbh,'box','off')
set(cbh,'LineWidth',LW)
set(gca,'LineWidth',LW)
set(cbh,'xtick',ceil(min(min(pktrampl(find(pktrampl~=0))))):ceil((floor(max(max(pktrampl)))-ceil(min(min(pktrampl(find(pktrampl~=0))))))/10):floor(max(max(pktrampl))))
set(cbh,'xticklabel',ceil(min(min(pktrampl(find(pktrampl~=0))))):ceil((floor(max(max(pktrampl)))-ceil(min(min(pktrampl(find(pktrampl~=0))))))/10):floor(max(max(pktrampl))))

set(gca, 'LooseInset', [0,0.2,0,0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Coefficient of Variation Figure%%%%%%%%%%%%%%
fh = figure;

set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
Pwidth = 0.45*Ssize(3);
Pheight = 1*Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

fh = pcolor(kgrid,Fgrid,CVpk);
colormap(ametrine);

cd=get(fh,'cdata');
cd(find(CVpk==0))=NaN;
set(fh,'cdata',cd)

%The color range is rescaled to match the nonzero data range
caxis([min(min(CVpk(find(CVpk~=0)))) max(max(CVpk))])
axis square;
xlabel({'Stiffness (mN\cdotm^{-1})' ''},'FontSize',FS,'FontName',FN);
ylabel('Force (pN)','FontSize',FS,'FontName',FN);
kticks = ((kgrid+circshift(kgrid,[1 -1]))/2)*10^1;
kticks(end) = [];
Fticks = ((Fgrid+circshift(Fgrid,[1 -1]))/2)*10^1;
Fticks(end) = [];
%Set the labels at the correct positions and round the label values
set(gca,'xtick',kticks,'ytick',Fticks,'FontSize',FS,'FontName',FN)
set(gca,'xticklabel',round(kticks),'ytick',round(Fticks),'FontSize',FS,'FontName',FN)

set(gca,'TickDir','Out')
set(gca, 'TickLength', [TL TL]);
set(gca,'box','off')
set(gca,'LineWidth',LW)
set(fh,'LineWidth',LW)

cbh = colorbar('location','Southoutside');
xlabel(cbh,'Amplitude (nm)','FontSize',FS,'FontName',FN);
set(gca,'FontSize',FS,'FontName',FN)
set(cbh, 'TickLength', [TL TL]);
set(cbh,'TickDir','Out')
set(cbh,'box','off')
set(cbh,'LineWidth',LW)
set(gca,'LineWidth',LW)
set(cbh,'xtick',ceil(min(min(CVpk(find(CVpk~=0))))):ceil((floor(max(max(CVpk)))-ceil(min(min(CVpk(find(CVpk~=0))))))/10):floor(max(max(CVpk))))
set(cbh,'xticklabel',ceil(min(min(CVpk(find(CVpk~=0))))):ceil((floor(max(max(CVpk)))-ceil(min(min(CVpk(find(CVpk~=0))))))/10):floor(max(max(CVpk))))

set(gca, 'LooseInset', [0,0.2,0,0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
