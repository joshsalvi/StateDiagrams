%Experimental data. All points files%%%%%%%%%%%%%%%%%
%load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-05.01/Ear 1/Cell 7/20140805-cell7.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-05.01/Ear 1/Cell 7/Tstartend2Hzmin.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/20130908-cell15-4-2d.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/4-Tstartend1sec1Hzmin.mat');

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
tmin = 0;%ms
tmax = tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)));

%%%%%%%%%%%%%%%%%%%%Steady-state%%%%%%%%%%%
%Check that the local mean and std have converged

%Remove the drift from each time trace
fdrift = 0.001;%kHz
fdrift = 0.0005;%kHz            % CHOICE
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

Countmin = 0;
Countmax = 4.5*10^-2;
Xmax = 20;%nm
Xmin = -Xmax;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Grid Loops%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Findex = 1:length(Fgrid)
for kindex = 1:length(kgrid)
%for Findex = 1:2
%for kindex = 1:2
        
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

%%%%%%%%%%%%%Displacement Histograms%%%%%%%%%%    
    
%wsdetrend should be independent of each tt length
wsdetrend = round(ws);
X = Xd(ssstartpt:ssendpt,Npulse,Ntrial)- smooth(Xd(ssstartpt:ssendpt,Npulse,Ntrial),wsdetrend,'moving');

%%%%%%%%Symmetry test%%%%%%%%%%%%%
%Xctr = median(X) sometimes => rejection symmetric bimodal distributions
Xctr = mean(X);

%Asai06
%Asymmetric if Xupper doesn't come from same distribution as Xlower
Xupper = X(X>=Xctr) - Xctr;
Xlower = Xctr - X(X<=Xctr);

%Asai06 doesn't symmetrize. Symmetrize => more likely to accept the null.
Xupper = [Xupper' -Xupper']';
Xlower = [Xlower' -Xlower']';

%Null is the distribution is symmetric so the halfs are the same
%CHOICE
length(num2str(X(10),'%10.100g')); %Number of significant digits in data
% CHOICE
alphasymm = 10^-2;
eps(alphasymm); %The difference between alphasymm and the next largest double-precision number
[hsymm,psymm,KSstat] = kstest2(Xupper,Xlower,alphasymm);
% CHOICE
if psymm <= alphasymm && KSstat >= 0.02
    hsymm = 1;
    symmtext = 'A';
else
    hsymm = 0;
    symmtext = 'S';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Kurtosis test%%%%%%%%%%%%%%%%%%%%%
% CHOICE
alphaK = 10^-3;
NK = length(X);
%Reference: Basic Statistics for Social Research (1997)
%Standard error in the skewness
SES = sqrt(6*NK*(NK-1)/((NK-2)*(NK+1)*(NK+3)));
%Standard error in the kurtosis = standard deviation of the kurtosis
SEK = 2*SES*sqrt((NK^2-1)/((NK-3)*(NK+5)));
Kstat = (kurtosis(X,0)-3)/SEK;
%One-sided test to find distributions with very small kurtosis values
%Kstat is approximately normally distributed for NK > 10^3 with
%mean 0 and standard deviation 1
pK = cdf('Normal',Kstat,0,1);
% CHOICE
if pK <= alphaK && Kstat <= -18.0
    Ktext = 'F';
    hk = 1;
else
    Ktext = 'T';
    hk = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Unimodality test%%%%%%%%%%
%CHOICE
alphauni = 0.01;
%Null is that there is no dip. Null is unimodality.
%Algorithm is sort X and look for a concave/convex change
%dip ~> 0.1 => not unimodal
%Nboot = 5*10^4;
% CHOICE
Nboot = 2*10^2; %Number of bootstraped distributions needed to find puni
[dip, puni, Xlow, Xup]=hartigansdipsigniftest(X,Nboot);
% CHOICE
if puni <= alphauni && dip >= 0.002
    huni = 1;
    unitext = 'M';
else
    huni = 0;
    unitext = 'U';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bin width using the Freedman/Diaconis rule
bw = 2*iqr(X)/length(X)^(1/3);
Nbins = ceil((max(X) - min(X))/bw);
BinSize = (max(X) - min(X))/Nbins;
[binFreq,binPos]=hist(X,Nbins);

CountTotal = sum(binFreq);

%Normalize the histogram
binFreq = binFreq./CountTotal;

sph = subplot(length(Fgrid),length(kgrid),kindex+(Findex-1)*length(Fgrid));
bar(binPos,binFreq,'b');
hold on
stem([Xlow Xup],[binFreq(find(abs(Xlow-binPos)==min(abs(Xlow-binPos)))) binFreq(find(abs(Xup-binPos)==min(abs(Xup-binPos))))],'r')
%Distribution must be either asymmetric or fat to indicate oscillations
if hk == 1
set(gca,'Color',[0.5 1 0.5]);%Light Green
Mod(Npulse,Ntrial) = 1;
end
if hsymm == 1
set(gca,'Color',[1 1 0]);%Yellow
Mod(Npulse,Ntrial) = 1;
end
if huni == 1
set(gca,'Color',[1 0.8 0.8]);%Pink
Mod(Npulse,Ntrial) = 1;
end
axis([Xmin,Xmax,Countmin,Countmax])
SK = 0;
if SK == 1
if kurtosis(X) < 2.5
set(gca,'Color','y');
end
text(0.5*Xmin,0.8*Countmax,{'Skewness' num2str(skewness(X),'%3.1f')},...
    'FontSize',18,'HorizontalAlignment','center','Color','c');
text(0.5*Xmax,0.8*Countmax,{'Kurtosis' num2str(kurtosis(X),'%3.1f')},...
    'FontSize',18,'HorizontalAlignment','center','Color','m');
end
text(0.65*Xmax,0.25*Countmax,{num2str(pK,'%3.1e') num2str(Kstat,'%3.1e') Ktext},...
    'FontSize',18,'HorizontalAlignment','center','Color',[0 0.5 0]);
text(0.65*Xmax,0.75*Countmax,{num2str(psymm,'%3.1e') num2str(KSstat,'%3.1e') symmtext},...
    'FontSize',18,'HorizontalAlignment','center','Color',[1 0.3 0]);
text(0.6*Xmin,0.75*Countmax,{num2str(puni,'%3.1e') num2str(dip,'%3.1e') unitext},...
    'FontSize',18,'HorizontalAlignment','center','Color',[0.5 0 0]);
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
%%%%%%%%%%%Save the modality%%%%%%%%%%%
%modfile = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-08-05.01/Ear 1/Cell 7/Modality2Hzmin.mat';
modfile = '/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/SSOverTime/4-Modality1sec1Hzmin.mat';
save(modfile, 'Mod');
display('saving...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
