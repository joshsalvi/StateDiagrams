load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-04.01/Ear 1/Cell 6/20140804-cell6.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-04.01/Ear 1/Cell 6/Tstartend.mat');


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

deltat = 0.1;%ms
Fs = 1/deltat;%kHz
tvec = 0:deltat:(length(Xd)-1)*deltat;
tmin = 1000;%ms
tmax = 3000;
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)));


Countmin = 0;
Countmax = 4.5*10^-2;
Xmax = 10;%nm
Xmin = -20;

fh = figure;
set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
set(gca, 'LooseInset', get(gca, 'TightInset'));
Pwidth = Ssize(3);
Pheight = Ssize(4);
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

Fgridrev = sort(Fgrid,'descend');

%Modality = 1 => limit cycle oscillations 
Mod(1:Np,1:Nt) = zeros(Np,Nt);

h = -1;
matlabtest = 0;


%Remove the drift from each time trace
fdrift = 0.001;%Hz
fdrift = 0.0005;%Hz
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Tstartend(1,Npulse,Ntrial) ~= Tstartend(2,Npulse,Ntrial)
    ssstartpt = find(abs(tvec-Tstartend(1,Npulse,Ntrial))==min(abs(tvec-Tstartend(1,Npulse,Ntrial))));
    ssendpt = find(abs(tvec-Tstartend(2,Npulse,Ntrial))==min(abs(tvec-Tstartend(2,Npulse,Ntrial))));
    else
    ssstartpt = 1;
    ssstartpt = round(length(Xd)/2);
    ssendpt = length(Xd);
    end

    
    ssstartpt = 1;
    ssendpt = length(Xd);
    
    
%wsdetrend should be independent of each tt length
wsdetrend = round(ws);
X = Xd(ssstartpt:ssendpt,Npulse,Ntrial)- smooth(Xd(ssstartpt:ssendpt,Npulse,Ntrial),wsdetrend,'moving');
 

%Mean std
Xmeanstd = std(X(ssstartpt:ssendpt));


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
amplgrid(Findex,kindex) = Xmeanstd;
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
