%Finds peak of raw spectrum
%Use zero padding to get better estimate of peak amplitude
%Need to run hist2 first if clear all is commented out
%clear all
%close all

load('/Users/joshsalvi/Documents/Lab/Lab/Data Analysis/State Spaces/20130821/20130821-cell11.mat');

%time.data(:,1,1) is the time vec. Don't know what the other elements of
%time are but they're not time values. 

sizeXd = size(Xd_pulse);
Np = sizeXd(2);
Nt = sizeXd(3);

%Three points are duplicates at (0,0) 
Fsort = sort(F_rand);
Fgrid = Fsort(diff(Fsort) ~= 0);
Fgrid(end+1) = max(F_rand);
ksort = sort(k_rand);
kgrid = ksort(diff(ksort) ~= 0);
kgrid(end+1) = max(k_rand);

freq = zeros(length(Fgrid),length(kgrid));
ampl = zeros(length(Fgrid),length(kgrid));
rms = zeros(length(Fgrid),length(kgrid));

Fcpt = zeros(Np,Nt);
kpt = zeros(Np,Nt);

freq_rand = zeros(size(F_rand));
ampl_rand = zeros(size(F_rand));
rms_rand = zeros(size(F_rand));

%Set pl = 1 to output plots
pl = 0;
for Npulse = 1:Np
    for Ntrial = 1:Nt

AfterRamp = find(time.data(:,1,1) >= ramp);
BeforeRamp = find(time.data(:,1,1) <= time.data(end,1,1)-ramp);
%Choice
AfterRamp(1) = round(0.5*length(time.data(:,1,1)));

XdR = Xd(AfterRamp(1):BeforeRamp(end),Npulse,Ntrial);
timeR = time.data(AfterRamp(1):BeforeRamp(end),1,1);

%detrended data
%Choice
XdRS = XdR - smooth(XdR,round(length(XdR)/10),'sgolay',1);

if pl == 1
figure
plot(timeR,XdRS,'r')
title(['Bundle Displacement: k = ' num2str(10^6*k_rand(Npt)) ' mN/m, Fc = ' num2str(10^12*F_rand(Npt)) ' pN']); xlabel('Time (ms)'); ylabel('X (nm)');
end

%Fourier transform
%Choice
NFFT = numel(timeR);
NFFT2 = (2^4)*2^nextpow2(numel(timeR));
f = (Fs/2)*linspace(0,1,NFFT/2+1);
f2 = (Fs/2)*linspace(0,1,NFFT2/2+1);
effective_delta_f = Fs/NFFT2;
actual_delta_f = Fs/numel(timeR);
%double sided
XdRFT = fft(XdR,NFFT)/numel(timeR);
XdRFT2 = fft(XdR,NFFT2)/numel(timeR);
%single sided
XdRFTss = 2*XdRFT(1:length(f));
XdRFTss(1) = XdRFTss(1)/sqrt(2); %0 freq. component increases by sqrt(2)
XdRFTss(end) = XdRFTss(end)/sqrt(2); %last freq. component increases by sqrt(2)

XdRFTss2 = 2*XdRFT2(1:length(f2));
XdRFTss2(1) = XdRFTss2(1)/sqrt(2); %0 freq. component increases by sqrt(2)
XdRFTss2(end) = XdRFTss2(end)/sqrt(2); %last freq. component increases by sqrt(2)

if pl == 1
figure
plot(f,abs(XdRFTss),'red')
hold on
plot(f2,abs(XdRFTss2),'k')
set(gca,'XScale','log')
set(gca,'YScale','log')
end

minf = 5*actual_delta_f;
XdRmodspec = XdRFTss.*(f>minf)';
fpt = find(XdRmodspec == max(XdRmodspec));
delta_f_win = actual_delta_f;
fpoint_min = find( abs(f-(f(fpt)-delta_f_win)) == min(abs(f-(f(fpt)-delta_f_win))) );
fpoint_max = find( abs(f-(f(fpt)+delta_f_win)) == min(abs(f-(f(fpt)+delta_f_win))) );
XdRFTmodspecSTD = std(abs(XdRmodspec(round(fpoint_min/1):1*fpoint_max)));
XdRFTmodspecM = mean(abs(XdRmodspec(round(fpoint_min/1):1*fpoint_max)));

%RMS amplitude of detrended data
XdRSm = mean(XdRS);
XdRSrms = sqrt(sum((XdRS-XdRSm).^2)/length(XdRS));

Npt = (Ntrial-1)*Np+Npulse;
%Count the peak if it's signifcantly larger than the background
%if (max(abs(XdRmodspec(fpoint_min:fpoint_max)))-XdRFTmodspecM > 0*XdRFTmodspecSTD)
%Count the peak only if ampl dist is multimodal. Must run hist2 first.
if modality_rand(Npt) == 1
%freq_rand(Npt) = f(find(abs(XdRFTss) == max(abs(XdRFTss(fpoint_min:fpoint_max)))));
%ampl_rand(Npt) = max(abs(XdRFTss(fpoint_min:fpoint_max)));
%freq(find(Fgrid == F_rand(Npt)),find(kgrid == k_rand(Npt))) = f(find(abs(XdRFTss) == max(abs(XdRFTss(fpoint_min:fpoint_max)))));
%ampl(find(Fgrid == F_rand(Npt)),find(kgrid == k_rand(Npt))) = max(abs(XdRFTss(fpoint_min:fpoint_max)));

fmin_diff = min(abs(f2-f(fpoint_min)));
f2point_min = find((f2 == f(fpoint_min)+fmin_diff) | (f2 == f(fpoint_min)-fmin_diff));
fmax_diff = min(abs(f2-f(fpoint_max)));
f2point_max = find((f2 == f(fpoint_max)+fmax_diff) | (f2 == f(fpoint_max)-fmax_diff));
freq_rand(Npt) = f2(find(abs(XdRFTss2) == max(abs(XdRFTss2(f2point_min:f2point_max)))));
ampl_rand(Npt) = max(abs(XdRFTss2(f2point_min:f2point_max)));
freq(find(Fgrid == F_rand(Npt)),find(kgrid == k_rand(Npt))) = f2(find(abs(XdRFTss2) == max(abs(XdRFTss2(f2point_min:f2point_max)))));
ampl(find(Fgrid == F_rand(Npt)),find(kgrid == k_rand(Npt))) = max(abs(XdRFTss2(f2point_min:f2point_max)));
end
rms_rand(Npt) = XdRSrms;
rms(find(Fgrid == F_rand(Npt)),find(kgrid == k_rand(Npt))) = XdRSrms;
if pl == 1
    title(['Bundle Displacement: k = ' num2str(10^6*k_rand(Npt)) ' mN/m, Fc = ' num2str(10^12*F_rand(Npt)) ' pN']); xlabel('Time (ms)'); ylabel('X (nm)');
    text(0.1,2*10^-2*ampl_rand(Npt),['Peak Frequency (Hz) =' num2str(freq_rand(Npt))],'HorizontalAlignment','left')
    text(0.1,10^-2*ampl_rand(Npt),['Amplitude (nm) =' num2str(ampl_rand(Npt))],'HorizontalAlignment','left')
    text(0.1,5*10^-3*ampl_rand(Npt),['RMS Amplitude (nm) =' num2str(rms_rand(Npt))],'HorizontalAlignment','left')
end

    end
end

%Create grid so that data squares are centered correctly.
dFgrid = diff(Fgrid);
dFgrid(end+1) = dFgrid(end);
Fgridold = Fgrid;
Fgrid = Fgrid-dFgrid/2;
Fgrid(end+1) = Fgrid(end) + dFgrid(end);
freqold = freq;
amplold = ampl;
freq(end+1,:)=zeros(1,length(kgrid));
ampl(end+1,:)=zeros(1,length(kgrid));
rms(end+1,:)=zeros(1,length(kgrid));

dkgrid = diff(kgrid);
dkgrid(end+1) = dkgrid(end);
kgridold = kgrid;
kgrid = kgrid-dkgrid/2;
kgrid(end+1) = kgrid(end) + dkgrid(end);
freq(:,end+1)=zeros(length(Fgrid),1);
ampl(:,end+1)=zeros(length(Fgrid),1);
rms(:,end+1)=zeros(length(Fgrid),1);

kgrid=min(kgrid):(max(kgrid)-min(kgrid))/(length(kgrid)-1):max(kgrid);


figure
pcolor(kgrid*10^6,Fgrid*10^12,freq)
axis square;
title('Frequency (Hz)'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');
caxis([min(min(freq)) max(max(freq))])
colorbar
%shading interp

figure
pcolor(kgrid*10^6,Fgrid*10^12,ampl)
axis square;
title('Amplitude (nm)'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');
caxis([min(min(ampl)) max(max(ampl))])
colorbar
%shading interp

figure
pcolor(kgrid*10^6,Fgrid*10^12,rms)
axis square;
title('RMS Amplitude (nm)'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');
caxis([min(min(ampl)) max(max(ampl))])
colorbar
%shading interp

figure
scatter(Fgridold*10^12,freqold(:,1),'filled')
axis([min(Fgridold*10^12-1) max(Fgridold*10^12) 0 max(freqold(:,1)+1)])

%Look for correlation between ampl and freq
figure
hold on
freqsamp = freqold;
freqsamp(find(freqsamp==0)) = [];
amplsamp = amplold;
amplsamp(find(amplsamp==0)) = [];
scatter(freqsamp,amplsamp,'filled')
title('Amplitude vs Frequency'); xlabel('Frequency (Hz)'); ylabel('Amplitude (nm)');
corrmat = corrcoef(freqsamp,amplsamp);
text(0.7*max(freqsamp),0.9*max(amplsamp),['Correlation = ' num2str(corrmat(1,2))],'HorizontalAlignment','left')

for i=1:length(Fgridold)
freqsamp = freqold(i,:);
freqsamp(find(freqsamp==0)) = [];
amplsamp = amplold(i,:);
amplsamp(find(amplsamp==0)) = [];
figure
scatter(freqsamp,amplsamp,'filled')
corrmat = corrcoef(freqsamp,amplsamp);
if numel(corrmat) > 1
text(0.7*max(freqsamp),0.9*max(amplsamp),['Correlation = ' num2str(corrmat(1,2))],'HorizontalAlignment','left')
end
title(['Amplitude vs Frequency: F = ' num2str(Fgridold(i)*10^12)]); xlabel('Frequency (Hz)'); ylabel('Amplitude (nm)');
end

for i=1:length(kgridold)
freqsamp = freqold(:,i);
freqsamp(find(freqsamp==0)) = [];
amplsamp = amplold(:,i);
amplsamp(find(amplsamp==0)) = [];
figure
scatter(freqsamp,amplsamp,'filled')
corrmat = corrcoef(freqsamp,amplsamp);
if numel(corrmat) > 1
text(0.7*max(freqsamp),0.9*max(amplsamp),['Correlation = ' num2str(corrmat(1,2))],'HorizontalAlignment','left')
end
title(['Amplitude vs Frequency: k = ' num2str(kgridold(i)*10^6)]); xlabel('Frequency (Hz)'); ylabel('Amplitude (nm)');
end

%figure
%scatter3(k_rand*10^6,F_rand*10^12,freq_rand,100*ones(size(freq_rand)),freq_rand,'filled')
%view([0 90])
%title('Frequency (Hz)'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');

%figure
%scatter3(k_rand*10^6,F_rand*10^12,ampl_rand,100*ones(size(ampl_rand)),ampl_rand,'filled')
%view([0 90])
%title('Amplitude (nm)'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');

%figure
%scatter3(k_rand*10^6,F_rand*10^12,rms_rand,100*ones(size(rms_rand)),rms_rand,'filled')
%view([0 90])
%title('RMS Amplitude (nm)'); xlabel('stiffness (mN\cdotm^{-1})'); ylabel('force (pN)');
