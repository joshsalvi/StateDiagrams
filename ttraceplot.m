function ttraceplot(kindex,Findex,tmanstart,tmanend)
%Experimental data. All points files%%%%%%%%%%%%%%%%%
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-04.01/Ear 1/Cell 3/20140804-cell3.mat');
%load('/Users/dmelody/Work/Ear/Hair Bundle Expts/Models/FishAnal/smallcell2.mat');

%There are discontinuities at the ends of the traces
%k < 425 uN/m => G < 0

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

deltat = 0.1;%ms
Fs = 1/deltat;%kHz
tvec = 0:deltat:(length(Xd)-1)*deltat;
tmin = 0;%ms
tmax = tvec(end);

Npt = find(k_rand == kgrid(kindex) & F_rand == Fgridrev(Findex));
if rem(Npt,Np) ~= 0;
    Npulse = rem(Npt,Np);
else
    Npulse = Np;
end  
Ntrial = (Npt - Npulse)/Np + 1;
Npt = (Ntrial-1)*Np+Npulse;

%Remove the mean from each time trace
X = Xd(:,Npulse,Ntrial)-mean(Xd(:,Npulse,Ntrial));

%%%%%%%%%%%Save the time limits%%%%%%%%%%%
TLsave = 0;
if TLsave == 1
timefile = '/Users/dmelody/Work/Ear/Hair Bundle Expts/Models/Noise/Tstartend.mat';
%If the file exists then load it
if exist(timefile,'file') == 2
load(timefile);
display('loading...');
else %If file doesn't exist then create variable
Tstartend(1:2,1:Np,1:Nt) = zeros(2,Np,Nt);
display('creating...');
end
Tstartend(1,Npulse,Ntrial) = tmanstart;
Tstartend(2,Npulse,Ntrial) = tmanend;
save(timefile, 'Tstartend');
display('saving...');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
set(fh,'PaperPositionMode','auto')
Ssize =  get(0,'ScreenSize');
set(gca, 'LooseInset', get(gca, 'TightInset'));
Pwidth = 0.5*Ssize(3);
Pheight = Pwidth;
set(fh, 'Position', [(Ssize(3)-Pwidth)/2 (Ssize(4)-Pheight)/2 Pwidth Pheight])

TitleFS = 18;
TextFS = 16;
asp_rat = [0.8*Ssize(3)/Ssize(4) 1 1];

leftfac = 0.035;
bottomfac = 0.17;
widthfac = 0.13;
heightfac = 0.6;

findss = 0;
%%%%%%%%%%%%%%%%%%%%Steady-state%%%%%%%%%%%
%Check that the local mean and std have converged

%Remove the drift from each time trace
%CHOICE
fdrift = 0.001;%Hz
%CHOICE
fsmooth = 3*fdrift;
%Must smooth for a time less than 1/fdrift to remove fdrift
%ws must be odd
ws = (1/fsmooth)/deltat;
if rem(floor(ws),2) == 1 %odd
    ws = floor(ws);
else %even
    ws = floor(ws)+1; 
end

%Find the running mean and std for a time greater than the longest period (1/fsmooth)
%CHOICE
tss = 1.5*(1/fsmooth);
ssws = round(tss/deltat);
%CHOICE
nstdss = 2;
%CHOICE
stdfrac = 0.1;%Get to 90% of std ss
%Discount end discontinuities
%CHOICE
ndev = 2;

    if findss == 1
    %Discount end discontinuities
    Xdmad = mad(Xd(:,Npulse,Ntrial),1);
    Xdindex = length(Xd);
    if sum(abs(diff(Xd(:,Npulse,Ntrial))) > ndev*Xdmad) >= 1
    Xdindex = find(abs(diff(Xd(:,Npulse,Ntrial))) > ndev*Xdmad);
    end
    if min(Xdindex) >= length(Xd) - ssws + 1
    ssendpt = min(Xdindex);
    end
    end
    
    %Find the ss manually
    if findss == 0
    ssstartpt = find(abs(tvec-tmanstart)==min(abs(tvec-tmanstart)));
    ssendpt = find(abs(tvec-tmanend)==min(abs(tvec-tmanend)));
    end
    
    %Xnodrift is defined from 1 up to ssendpt
    %Moving average - first zero in FT at 1/wsdetrend
    %Moving is best for retaining sudden transitions
    %sgolay smooths out transitions and obscures oscillations more
    Xnodrift = Xd(1:ssendpt,Npulse,Ntrial) - smooth(Xd(1:ssendpt,Npulse,Ntrial),ws,'moving');
    
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
  
    if findss == 1
    ssstartpt = 1;
    for i = 1:length(Xmean)-length(Xsspts)
        %SS starts once the mean absolute deviation of the mean and std are within 
        %some % of their ss values
        if mean(abs(Xmean(length(Xmean)-length(Xsspts)+1-i:length(Xmean)-i)-Xmeanss))>nstdss*Xstdss || mean(abs(Xstd(length(Xstd)-length(Xsspts)+1-i:length(Xstd)-i)-Xstdss))>stdfrac*Xstdss

           %Previous window start is initial ss start estimate
           ssstartpt = length(Xmean)-length(Xsspts)+1-i+1;
           %Find the earliest point within the previous window that has mean & std ~ ss
           %May miss early part of ss
           ssstartpt = max([min(find(abs(Xmean(length(Xmean)-length(Xsspts)+1-i:length(Xmean)-i)-Xmeanss)<nstdss*Xstdss)) min(find((abs(Xstd(length(Xstd)-length(Xsspts)+1-i+1:length(Xstd)-i+1)-Xstdss)<stdfrac*Xstdss)))])+ssstartpt;
            break
        end
    end
    end
    
    %Mean local std
    Xmeanstd = mean(Xstd(ssstartpt:end));

%%%%%%%%%%%SS Time Trace Plot%%%%%%%%
Xmax = 1.1*max(X);
Xmin = 1.1*min(X);

%%%%%%%%%%No drift plot%%%%%%%%%%%%%
%Xnodrift is only defined up to ssendpt
Xmax = max(Xnodrift(ssstartpt:end));
Xmin = min(Xnodrift(ssstartpt:end));
plot(tvec(ssstartpt:ssendpt)-tvec(ssstartpt),Xnodrift(ssstartpt:end),'g')
axis([tvec(ssstartpt)-tvec(ssstartpt),tvec(ssendpt)-tvec(ssstartpt),Xmin,Xmax])
pbaspect(asp_rat)
set(gca,'LineWidth',2)
set(gca,'FontSize',36)
text((tvec(ssendpt)+tvec(ssstartpt))/2,0,{'Local RMS' [num2str(Xmeanstd,'%4.1f') ' nm']},...
    'FontSize',TextFS,'HorizontalAlignment','center');

ylabel('Displacement (nm)','FontSize',36);
set(gca,'ytickMode', 'auto')
xlabel('Time (ms)','FontSize',36);
set(gca,'xtickMode', 'auto')
