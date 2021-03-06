%Experimental data. All points files%%%%%%%%%%%%%%%%%
%clear;close all;clc;
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Extracted Data-good3.mat')
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Tstartend_dipstat.mat');
figure;
fishfigs = 1;
if fishfigs == 1
%load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/Gentamicin/2014-08-05.01/Ear 1/Cell 8/Modality-foranalysis.mat');
modfile = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Modality2sec2Hzmin-good-dipstat-filtered.mat';
load(modfile);
load('/Users/joshsalvi/GitHub/StateDiagrams/customcolormaps-redblue.mat');
end

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
amplgrid = zeros(length(Fgrid),length(kgrid));

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

%Remove the mean
X = Xd(ssstartpt:ssendpt,Npulse,Ntrial)-mean(Xd(ssstartpt:ssendpt,Npulse,Ntrial));
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
% Eliminate 60,120,180 � 1 Hz
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
    %freqrange= (fpsd <= fmax & fpsd >= fmin & Xpsd > circshift(Xpsd,[1 1]) & Xpsd > circshift(Xpsd,[-1 1]));
    freqrange= (fpsd <= fmax & fpsd >= fmin & Xpsd);
end

Xpsdmaxima = Xpsd.*(freqrange);

Xpsdpeak = max(Xpsdmaxima);

XpsdmaximaMedian = median(Xpsdmaxima(find(Xpsdmaxima >= 0.1*Xpsdpeak)));

%The main peak of bundle oscillations. Can be bigger than the local std.
XFTpeak = (sqrt(fscale.*Xpsdpeak.*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
fpsdpeak = fpsd(find(Xpsd == Xpsdpeak));
if numel(fpsdpeak) > 1
    fpsdpeak = fpsdpeak(1);
    XFTpeak = XFTpeak(1);
elseif isempty(fpsdpeak) == 1
    fpsdpeak = fpsd(find(Xpsd==max(Xpsd)));
    Xpsdpeak = Xpsd(find(Xpsd==max(Xpsd)));
    XFTpeak = (sqrt(fscale.*Xpsdpeak.*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
end

sph = subplot(length(Fgrid),length(kgrid),kindex+(Findex-1)*length(Fgrid));
plot(fpsd,Xpsd,'m')
hold on
plot(fpsdpeak,Xpsdpeak,'k.','markersize',10);
if Tstartend(1,Npulse,Ntrial) ~= Tstartend(2,Npulse,Ntrial)
text(0.5*fmax,0.5*Xpsdmaxlim,{[num2str(Tstartend(1,Npulse,Ntrial)/1000,'%3.1f') '-' num2str(Tstartend(2,Npulse,Ntrial)/1000,'%3.1f') ' s']},...
    'FontSize',12,'HorizontalAlignment','center');
end
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([fmin,fmax,Xpsdminlim,Xpsdmaxlim])
%}
%text(fmax/2,0.8*Xpsdmax,{num2str([kgrid(kindex)*10^1,Fgridrev(Findex)*10^1])},...
%    'FontSize',12,'HorizontalAlignment','center');
text(0.03,0.1*Xpsdmaxlim,{[num2str(XFTpeak,'%3.1f') ' nm   ' num2str(fpsdpeak*fscale,'%3.1f') ' Hz']},...
    'FontSize',12,'HorizontalAlignment','center');
text(0.03,0.01*Xpsdmaxlim,{num2str(Xpsdpeak/XpsdmaximaMedian,'%3.1f')},...
    'FontSize',12,'HorizontalAlignment','center');

%If there are oscillations, record their frequency and amplitude
if Mod(Npulse,Ntrial) == 1
freqgrid(Findex,kindex) = fpsdpeak;%Convert to Hz
amplgrid(Findex,kindex) = XFTpeak;%nm
set(gca,'Color',[1 1 0]);%Yellow
end

spp = get(sph, 'pos');
set(sph, 'Position', [spp(1) spp(2) 1.3*spp(3) 1.4*spp(4)]);
    
if kgrid(kindex) == min(kgrid)
ylabel(num2str(Fgridrev(Findex)*10^1),'FontSize',12);
set(gca,'ytickMode', 'auto')
else
set(gca,'YTickLabel',[])
end
if Fgridrev(Findex) == min(Fgrid)
xlabel(num2str(kgrid(kindex)*10^1),'FontSize',12);
set(gca,'xtickMode', 'auto')
else
set(gca,'XTickLabel',[])
end
%}
end
end

if fishfigs == 0
return
end

%%%%%%%%%%%%%%%%%%%%Fish Figures%%%%%%%%%%%%%%

%Flip grid along the force direction
freqgrid = flipdim(freqgrid,1);
amplgrid = flipdim(amplgrid,1);

%%%%%Frequency and amplitude correlations%%%%%%%%%
amplvec = reshape(amplgrid,[length(amplgrid(:,1))*length(amplgrid(1,:)) 1]);
%Find the points that are within the fish
fishvecpts = find(amplvec ~=0);
amplvec = amplvec(fishvecpts);
freqvec = reshape(freqgrid,[length(freqgrid(:,1))*length(freqgrid(1,:)) 1]);
freqvec = freqvec(fishvecpts);

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

save('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/FreqAmplcorrelations2sec2Hzmin-good1.mat','rhofreqampl','prhofreqampl','rhofreqk','prhofreqk','rhofreqF','prhofreqF','rhoamplk','prhoamplk','rhoamplF','prhoamplF')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create grid so that data squares are centered correctly.
dFgrid = diff(Fgrid);
dFgrid(end+1) = dFgrid(end);
Fgridold = Fgrid;
Fgrid = Fgrid-dFgrid/2;
Fgrid(end+1) = Fgrid(end) + dFgrid(end);
freqgrid(end+1,:)=zeros(1,length(kgrid));
amplgrid(end+1,:)=zeros(1,length(kgrid));

dkgrid = diff(kgrid);
dkgrid(end+1) = dkgrid(end);
kgridold = kgrid;
kgrid = kgrid-dkgrid/2;
kgrid(end+1) = kgrid(end) + dkgrid(end);
freqgrid(:,end+1)=zeros(length(Fgrid),1);
amplgrid(:,end+1)=zeros(length(Fgrid),1);

kgrid=min(kgrid):(max(kgrid)-min(kgrid))/(length(kgrid)-1):max(kgrid);

FS = 24; %Font Size
FN = 'Arial'; %Font Name
LW = 2; %LineWidth
TL = 0.015; %TickLength

%%%%%%%%%%%%%%%%%%%%Frequency Figure%%%%%%%%%%%%%%
fh = figure;

set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
Pwidth = 0.45*Ssize(3);
Pheight = 1*Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

fh = pcolor(kgrid*10^1,Fgrid*10^1,freqgrid);
colormap(blue1);

cd=get(fh,'cdata');
cd(find(freqgrid==0))=NaN;
set(fh,'cdata',cd)

%The color range is rescaled to match the nonzero data range
caxis([min(min(freqgrid(find(freqgrid~=0)))) max(max(freqgrid))])
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
xlabel(cbh,'Frequency (Hz)','FontSize',FS,'FontName',FN);
set(gca,'FontSize',FS,'FontName',FN)
set(cbh, 'TickLength', [TL TL]);
set(cbh,'TickDir','Out')
set(cbh,'box','off')
set(cbh,'LineWidth',LW)
set(gca,'LineWidth',LW)
set(cbh,'xtick',ceil(min(min(freqgrid(find(freqgrid~=0))))):ceil((floor(max(max(freqgrid)))-ceil(min(min(freqgrid(find(freqgrid~=0))))))/10):floor(max(max(freqgrid))))
set(cbh,'xticklabel',ceil(min(min(freqgrid(find(freqgrid~=0))))):ceil((floor(max(max(freqgrid)))-ceil(min(min(freqgrid(find(freqgrid~=0))))))/10):floor(max(max(freqgrid))))

set(gca, 'LooseInset', [0,0.2,0,0]);

%saveas(gcf,'FreqFish','epsc')

%shading interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Amplitude Figure%%%%%%%%%%%%%%
fh = figure;

set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
Pwidth = 0.45*Ssize(3);
Pheight = 1*Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

fh = pcolor(kgrid*10^1,Fgrid*10^1,amplgrid);
colormap(red1);

cd=get(fh,'cdata');
cd(find(amplgrid==0))=NaN;
set(fh,'cdata',cd)

%The color range is rescaled to match the nonzero data range
caxis([min(min(amplgrid(find(amplgrid~=0)))) max(max(amplgrid))])
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
set(cbh,'xtick',ceil(min(min(amplgrid(find(amplgrid~=0))))):ceil((floor(max(max(amplgrid)))-ceil(min(min(amplgrid(find(amplgrid~=0))))))/10):floor(max(max(amplgrid))))
set(cbh,'xticklabel',ceil(min(min(amplgrid(find(amplgrid~=0))))):ceil((floor(max(max(amplgrid)))-ceil(min(min(amplgrid(find(amplgrid~=0))))))/10):floor(max(max(amplgrid))))

set(gca, 'LooseInset', [0,0.2,0,0]);

%saveas(gcf,'AmplFish','epsc')
%shading interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
