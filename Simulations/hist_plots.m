%Experimental data. All points files%%%%%%%%%%%%%%%%%
load('Tstartend1.0Noise.mat');
xfishtable = readtable('xfish1.0Noise.dat','Delimiter','\t');

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

%Really the probability limits
Countmin = 0;
Countmax = 1*10^-1;
Xmax = 4;
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
        
Npt = find(k_rand == kgrid(kindex) & F_rand == Fgridrev(Findex));
if rem(Npt,Np) ~= 0;
    Npulse = rem(Npt,Np);
else
    Npulse = Np;
end  
Ntrial = (Npt - Npulse)/Np + 1;
Npt = (Ntrial-1)*Np+Npulse;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ssstartpt = 1;
    ssendpt = length(Xd);

%%%%%%%%%%%%%Displacement Histograms%%%%%%%%%%    
    
%No need to detrend
X = Xd(ssstartpt:ssendpt,Npulse,Ntrial);
Kval = k_rand(Npulse,Ntrial);
Fcval = F_rand(Npulse,Ntrial);

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
alphasymm = 10^-3;
eps(alphasymm); %The difference between alphasymm and the next largest double-precision number
%KSstat is the max absolute difference btw the cdfs
[hsymm,psymm,KSstat] = kstest2(Xupper,Xlower,alphasymm);
if psymm <= alphasymm && KSstat >= 8*10^-2
    hsymm = 1;
    symmtext = 'A';
else
    hsymm = 0;
    symmtext = 'S';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Kurtosis test%%%%%%%%%%%%%%%%%%%%%
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
if pK <= alphaK && Kstat <= -60.0
    Ktext = 'F';
    hk = 1;
else
    Ktext = 'T';
    hk = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Unimodality test%%%%%%%%%%
%CHOICE
alphauni = 0.05;
%Null is that there is no dip. Null is unimodality.
%Algorithm is sort X and look for a concave/convex change
%dip ~> 0.1 => not unimodal
%Nboot = 5*10^4;
%Nboot = 1*10^2;
Nboot = 1*10^1; %Number of bootstraped distributions needed to find puni
[dip, puni, Xlow, Xup]=hartigansdipsigniftest(X,Nboot);
if puni <= alphauni && dip >= 1*10^-3
    huni = 1;
    unitext = 'M';
else
    huni = 0;
    unitext = 'U';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bin width using the Freedman/Diaconis rule
bw = 2*iqr(X)/length(X)^(1/3);
%Nbins = round((max(X) - min(X))/bw);
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

%K and Fc labels to check labeling
%text(0.6*Xmin,0.3*Countmax,{num2str(Kval,'%1.2f') num2str(Fcval,'%1.2f')},...
%    'FontSize',18,'HorizontalAlignment','center','Color','k');

%text(0.6*Xmin,0.3*Countmax,{num2str(convergeddist,'%3.1e') num2str(NbinsFrac,'%3.1e')},...
%    'FontSize',18,'HorizontalAlignment','center','Color','k');

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
%%%%%%%%%%%Save the modality%%%%%%%%%%%
modfile = 'Modality.mat';
save(modfile, 'Mod');
display('saving...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%