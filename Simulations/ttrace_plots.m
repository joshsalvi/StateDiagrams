%Scale the fish diagrams

xfishtable = readtable('xfish1.0Noise.dat','Delimiter','\t');
load('Modality1.0Noise.mat');
fishfigs = 1;
if fishfigs == 1
load('customcolormaps-redblue.mat');
end
load('Tstartend1.0Noise.mat');

time = xfishtable{3:end,{'time'}};

%Get the number of K and Fc points
Kno = xfishtable{1,{'time'}};
Fcno = xfishtable{2,{'time'}};

k_rand = zeros(Kno,Fcno);
F_rand = zeros(Kno,Fcno);
Xd = zeros(length(time),Kno,Fcno);
%File indices start at 0 but indices start at 1 in matlab
for Fcindex = 1:Fcno
    for Kindex = 1:Kno
k_rand(Kindex,Fcindex) = xfishtable{1,{['OP' num2str(Kindex-1) num2str(Fcindex-1)]}};
F_rand(Kindex,Fcindex) = xfishtable{2,{['OP' num2str(Kindex-1) num2str(Fcindex-1)]}};
Xd(:,Kindex,Fcindex) = xfishtable{3:end,{['OP' num2str(Kindex-1) num2str(Fcindex-1)]}};
    end
end

%Operating points in ascending order
Fsort = sort(F_rand(1,:));
Fgrid = Fsort(diff(Fsort) ~= 0);
Fgrid(end+1) = max(F_rand(1,:));
ksort = sort(k_rand(:,1));
kgrid = ksort(diff(ksort) ~= 0)';
kgrid(end+1) = max(k_rand(:,1));

sizeXd = size(Xd);
Np = sizeXd(2);
Nt = sizeXd(3);

deltat = time(2) - time(1);
Fs = 1/deltat;
tvec = time;
tmin = tvec(1);
tmax = tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)));

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
    
    ssstartpt = 1;
    ssendpt = length(Xd);
    
X = Xd(:,Npulse,Ntrial);

%Mean std
Xmeanstd = std(X(ssstartpt:ssendpt));

Kval = k_rand(Npulse,Ntrial);
Fcval = F_rand(Npulse,Ntrial);

sph = subplot(length(Fgrid),length(kgrid),kindex+(Findex-1)*length(Fgrid));
plot(tvec,X,'m')
hold on
plot(tvec(ssstartpt:ssendpt),X(ssstartpt:ssendpt),'g')
if Tstartend(1,Npulse,Ntrial) ~= Tstartend(2,Npulse,Ntrial)
text(0.5*tmax,0.8*Xmax,{[num2str(Tstartend(1,Npulse,Ntrial)/1000,'%3.1f') '-' num2str(Tstartend(2,Npulse,Ntrial)/1000,'%3.1f') ' s']},...
    'FontSize',12,'HorizontalAlignment','center');
end
axis([tmin,tmax,Xmin,Xmax])
text(tmax/2,0.6*Xmin,{'Local RMS' [num2str(Xmeanstd,'%4.1f') 'nm']},...
    'FontSize',12,'HorizontalAlignment','center');

%text(tmax/2,0.8*Xmax,{num2str([kgrid(kindex)*Kscale,Fgridrev(Findex)*Fcscale])},...
%    'FontSize',12,'HorizontalAlignment','center');

%K and Fc labels to check labeling
%text(0.1*tmax,0.5*Xmin,{num2str(Kval,'%1.2f') num2str(Fcval,'%1.2f')},...
%    'FontSize',18,'HorizontalAlignment','center','Color','k');


%If there are oscillations, record the rms amplitude
if Mod(Npulse,Ntrial) == 1
amplgrid(Findex,kindex) = Xmeanstd*amplscale;
set(gca,'Color',[1 1 0]);%Yellow
end

spp = get(sph, 'pos');
set(sph, 'Position', [spp(1) spp(2) 1.3*spp(3) 1.4*spp(4)]);

if kgrid(kindex) == min(kgrid)
ylabel(num2str(Fgridrev(Findex)),'FontSize',12);
set(gca,'ytickMode', 'auto')
else
set(gca,'YTickLabel',[])
end
if Fgridrev(Findex) == min(Fgrid)
xlabel(num2str(kgrid(kindex)),'FontSize',12);
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
save('RMScorrelations.mat','rhoRMSk','prhoRMSk','rhoRMSF','prhoRMSF')
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

fh = pcolor(kgrid*Kscale,Fgrid*Fcscale,amplgrid);

hold on
a = 3.5;
b = 0.5;
tau = 10;
kloopmax = a-1/tau;
kloop = min(kgrid):0.01:kloopmax;
hold on
plot(Kscale*kloop,real(Fcscale.*sqrt(((a-kloop-1./tau)./27).*(((2.*a+1./tau)*(1-b)-(2+b).*kloop)./(1-b)).^2)),'g','LineWidth',4)
plot(Kscale*kloop,-real(Fcscale.*sqrt(((a-kloop-1./tau)./27).*(((2.*a+1./tau)*(1-b)-(2+b).*kloop)./(1-b)).^2)),'g','LineWidth',4)

ktailmax = a*(1-b);
ktail = min(kgrid):0.01:ktailmax;
plot(Kscale*ktail,real(Fcscale.*sqrt((4/27).*((a.*(1-b)-ktail)./(1-b)).^3)),'g','LineWidth',4)
plot(Kscale*ktail,-real(Fcscale.*sqrt((4/27).*((a.*(1-b)-ktail)./(1-b)).^3)),'g','LineWidth',4)

colormap(red1);

cd=get(fh,'cdata');
cd(find(amplgrid==0))=NaN;
set(fh,'cdata',cd)

%The color range is rescaled to match the nonzero data range
caxis([min(min(amplgrid(find(amplgrid~=0)))) max(max(amplgrid))])
axis square;
xlabel({'Stiffness' ''},'FontSize',FS,'FontName',FN);
ylabel('Force','FontSize',FS,'FontName',FN);
kticks = ((kgrid+circshift(kgrid,[1 -1]))/2)*Kscale;
kticks(end) = [];
Fticks = ((Fgrid+circshift(Fgrid,[1 -1]))/2)*Fcscale;
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
xlabel(cbh,'Amplitude','FontSize',FS,'FontName',FN);
set(gca,'FontSize',FS,'FontName',FN)
set(cbh, 'TickLength', [TL TL]);
set(cbh,'TickDir','Out')
set(cbh,'box','off')
set(cbh,'LineWidth',LW)
set(gca,'LineWidth',LW)
set(cbh,'xtick',ceil(min(min(amplgrid(find(amplgrid~=0))))):ceil((floor(max(max(amplgrid)))-ceil(min(min(amplgrid(find(amplgrid~=0))))))/10):floor(max(max(amplgrid))))
set(cbh,'xticklabel',ceil(min(min(amplgrid(find(amplgrid~=0))))):ceil((floor(max(max(amplgrid)))-ceil(min(min(amplgrid(find(amplgrid~=0))))))/10):floor(max(max(amplgrid))))

set(gca, 'LooseInset', [0,0.2,0,0]);

saveas(gcf,'RMSFish','epsc')
%shading interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
