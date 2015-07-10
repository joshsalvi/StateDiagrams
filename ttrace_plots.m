clear all; close all;
%%%%%%%

load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Extracted Data.mat')
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Tstartend.mat');
fishfigs = 1;
if fishfigs == 1
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/Modality2sec2Hzmin.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/Gentamicin/2014-08-05.01/Ear 1/Cell 7/Modality2sec2Hzmin.mat');
load('/Users/joshsalvi/GitHub/StateDiagrams/customcolormaps-redblue.mat');
end
%}
%{
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/20130908-cell15-3-2d.mat');
fishfigs = 1;
if fishfigs == 1
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/3-Modality2sec1Hzmin.mat');
load('/Users/joshsalvi/GitHub/StateDiagrams/customcolormaps-redblue.mat');
end
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/3-Tstartend2sec1Hzmin.mat');
%}
%%%%%%%

%Operating points in ascending order
Fsort = sort(F_rand);
Fgrid = Fsort(diff(Fsort) ~= 0);
Fgrid(end+1) = max(F_rand);
ksort = sort(k_rand);
kgrid = ksort(diff(ksort) ~= 0);
kgrid(end+1) = max(k_rand);

sizeXd = size(Xd);
Np = sizeXd(2);
if size(sizeXd)==3
    Nt = sizeXd(3);
else
    Nt = 1;
end

deltat = 0.1;%ms
Fs = 1/deltat;%kHz
tvec = 0:deltat:(length(Xd)-1)*deltat;
tmin = 0;%ms
tmax = tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)));

%%%%%%%%%%%%%%%%%%%%Steady-state%%%%%%%%%%%
%Check that the local mean and std have converged

%Remove the drift from each time trace
fdrift = 0.002;%Hz
fsmooth = 3*fdrift;
%ws must be odd
ws = (1/fsmooth)/deltat;
if rem(floor(ws),2) == 1 %odd
    ws = floor(ws);
else %even
    ws = floor(ws)+1; 
end
%Find the running mean and std for a time greater than the longest period (1/fdrift)
tss = 1.5*(1/fsmooth);
ssws = round(tss/deltat);
nstdss = 2;
stdfrac = 0.1;%Get to 90% of std ss
%Discount end discontinuities
ndev = 2;

%Remove the mean from each time trace
Xd = bsxfun(@minus, Xd, mean(Xd(minindex:maxindex,:,:)));

Xmax = 2*max(max(std(Xd(minindex:maxindex,:,:))));
Xmin = -Xmax;

fh = figure;
set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
set(gca, 'LooseInset', get(gca, 'TightInset'));
Pwidth = Ssize(3);
Pheight = Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

Fgridrev = sort(Fgrid,'descend');

%Amplitude grid 
amplgrid = zeros(length(Fgrid),length(kgrid));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Grid Loops%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    if Tstartend(1,Npulse,Ntrial) ~= Tstartend(2,Npulse,Ntrial)
    ssstartpt = find(abs(tvec-Tstartend(1,Npulse,Ntrial))==min(abs(tvec-Tstartend(1,Npulse,Ntrial))));
    ssendpt = find(abs(tvec-Tstartend(2,Npulse,Ntrial))==min(abs(tvec-Tstartend(2,Npulse,Ntrial))));
    else
    ssstartpt = round(length(Xd)/2);
    ssendpt = length(Xd);
    end
    
    %%%%%%%%Steady State%%%%%%%%%%%%
    Xnodrift = Xd(1:ssendpt,Npulse,Ntrial) - smooth(Xd(1:ssendpt,Npulse,Ntrial),ws,'sgolay',1);
    tvec = 0:deltat:(length(Xd)-1)*deltat;
    
    %Xmean and Xstd are ssws-1 points shorter than Xnodrift 
    Xmean = zeros(1,length(Xnodrift)-ssws+1);
    Xstd = zeros(1,length(Xnodrift)-ssws+1);
    tvecss = zeros(1,length(Xnodrift)-ssws+1);
    for i = 1:length(Xnodrift)-ssws+1
        Xmean(i) = mean(Xd(i:ssws+i-1,Npulse,Ntrial));
        Xstd(i) = std(Xnodrift(i:ssws+i-1));
        tvecss(i) = (tvec(i)+tvec(ssws+i-1))/2.0;
    end

    %Define ss as a fixed time from the end of Xmean
    Xsspts = find(abs(tvecss-(tvecss(end)-tss))==min(abs(tvecss-(tvecss(end)-tss)))):length(Xmean);

    Xmeanss = mean(Xmean(Xsspts));
    Xstdss = mean(Xstd(Xsspts));
    
    %Mean local std
    Xmeanstd = mean(Xstd(ssstartpt:end));
    
X = Xd(:,Npulse,Ntrial);

sph = subplot(length(Fgrid),length(kgrid),kindex+(Findex-1)*length(Fgrid));
plot(tvec,X,'m')
hold on
plot(tvec(ssstartpt:ssendpt),X(ssstartpt:ssendpt),'g')
if Tstartend(1,Npulse,Ntrial) ~= Tstartend(2,Npulse,Ntrial)
text(0.5*tmax,0.8*Xmax,{[num2str(Tstartend(1,Npulse,Ntrial)/1000,'%3.1f') '-' num2str(Tstartend(2,Npulse,Ntrial)/1000,'%3.1f') ' s']},...
    'FontSize',12,'HorizontalAlignment','center');
end
axis([tmin,tmax,-50,50])
text(tmax/2,0.6*Xmin,{'Local RMS' [num2str(Xmeanstd,'%4.1f') 'nm']},...
    'FontSize',12,'HorizontalAlignment','center');

%text(tmax/2,0.8*Xmax,{num2str([kgrid(kindex)*10^6,Fgridrev(Findex)*10^12])},...
%    'FontSize',12,'HorizontalAlignment','center');

%If there are oscillations, record the rms amplitude
if Mod(Npulse,Ntrial) == 1
amplgrid(Findex,kindex) = Xmeanstd;%nm
set(gca,'Color',[1 1 0]);%Yellow
end

spp = get(sph, 'pos');
set(sph, 'Position', [spp(1) spp(2) 1.3*spp(3) 1.4*spp(4)]);

if kgrid(kindex) == min(kgrid)
ylabel(num2str(Fgridrev(Findex)*10^12),'FontSize',12);
set(gca,'ytickMode', 'auto')
else
set(gca,'YTickLabel',[])
end
if Fgridrev(Findex) == min(Fgrid)
xlabel(num2str(kgrid(kindex)*10^6),'FontSize',12);
set(gca,'xtickMode', 'auto')
else
set(gca,'XTickLabel',[])
end

end
end

if fishfigs == 0
return
end

%%%%%%%%%%%%%%%%%%%%Fish Figures%%%%%%%%%%%%%%

%Flip grid along the force direction
amplgrid = flipdim(amplgrid,1);

%%%%%%%%%%RMS correlations%%%%%%%%%%%%%%%%%%%%%%%%
amplvec = reshape(amplgrid,[length(amplgrid(:,1))*length(amplgrid(1,:)) 1]);
%Find the points that are within the fish
fishvecpts = find(amplvec ~=0);
amplvec = amplvec(fishvecpts);

Fgrid2 = Fgrid'*ones(1,length(Fgrid));
kgrid2 = ones(length(kgrid),1)*kgrid;
Fvec = reshape(Fgrid2,[length(Fgrid2(:,1))*length(Fgrid2(1,:)) 1]);
Fvec = Fvec(fishvecpts);
kvec = reshape(kgrid2,[length(kgrid2(:,1))*length(kgrid2(1,:)) 1]);
kvec = kvec(fishvecpts);
[rhoRMSk,prhoRMSk]=corr(amplvec,kvec,'type','Spearman');
[rhoRMSF,prhoRMSF]=corr(amplvec,Fvec,'type','Spearman');
%save('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-05.01/Ear 1/Cell 8/RMScorrelations2Hzmin.mat','rhoRMSk','prhoRMSk','rhoRMSF','prhoRMSF')
save('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-07-07.01/Ear 1/Cell 1/RMScorrelations2sec2Hzmin.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create grid so that data squares are centered correctly.
dFgrid = diff(Fgrid);
dFgrid(end+1) = dFgrid(end);
Fgridold = Fgrid;
Fgrid = Fgrid-dFgrid/2;
Fgrid(end+1) = Fgrid(end) + dFgrid(end);
amplgrid(end+1,:)=zeros(1,length(kgrid));

dkgrid = diff(kgrid);
dkgrid(end+1) = dkgrid(end);
kgridold = kgrid;
kgrid = kgrid-dkgrid/2;
kgrid(end+1) = kgrid(end) + dkgrid(end);
amplgrid(:,end+1)=zeros(length(Fgrid),1);

kgrid=min(kgrid):(max(kgrid)-min(kgrid))/(length(kgrid)-1):max(kgrid);

FS = 24; %Font Size
FN = 'Arial'; %Font Name
LW = 2; %LineWidth
TL = 0.015; %TickLength

%%%%%%%%%%%%%%%%%%%%Amplitude Figure%%%%%%%%%%%%%%
fh = figure;

set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
Pwidth = 0.45*Ssize(3);
Pheight = 1*Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

fh = pcolor(kgrid*10^6,Fgrid*10^12,amplgrid);
colormap(red1);

cd=get(fh,'cdata');
cd(find(amplgrid==0))=NaN;
set(fh,'cdata',cd)

%The color range is rescaled to match the nonzero data range
caxis([min(min(amplgrid(find(amplgrid~=0)))) max(max(amplgrid))])
axis square;
xlabel({'Stiffness (mN\cdotm^{-1})' ''},'FontSize',FS,'FontName',FN);
ylabel('Force (pN)','FontSize',FS,'FontName',FN);
kticks = ((kgrid+circshift(kgrid,[1 -1]))/2)*10^6;
kticks(end) = [];
Fticks = ((Fgrid+circshift(Fgrid,[1 -1]))/2)*10^12;
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

%saveas(gcf,'RMSFish','epsc')
%shading interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
