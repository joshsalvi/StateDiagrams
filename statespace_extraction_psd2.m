%% Import and Extract Data

 
clear all; close all;
 
iter_loc = 110; %logdata.data(i,iter_loc) : number of state points in the file i
alpha_loc=69;
ramp_loc1=1; %ramp_loc2=37;
path = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/';
date = '2014-08-04';
session = '.01';
ear = '1';
cell = '8';
user = 'JS';
sessions = ['152628 2';
'152826 3';
'153018 4';
'153158 5';
'153311 6'];
statespace = '1';
k_hb = 425e-6;
 
spacefile = sprintf('%s%s%s%s %s%s %s%s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/State space',statespace,'.txt');
spacedata=importdata(spacefile,'\t',3);
logfile = sprintf('%s%s%s%s %s%s %s%s%s %s%s %s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',user,date,session,ear,cell,'.log');
logdata = importdata(logfile);
 
spacedata.data=spacedata.data; %For some reason in some files, transpose it.
 
a = length(sessions(:,1));
 
for i = 1:a
file(i,:) = sprintf('%s%s%s%s %s%s %s%s%s %s%s %s %s %s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',user,date,session,ear,cell,'MS',sessions(i,:),'.txt');
end
 
for i = 1:a
data=importdata(file(i,:));
 
if i<a
    data2.data(:,:,i)=data.data;
else
    if length(data.data(:,1,1))<length(data2.data(:,1,1))
        data.data(length(data2.data(:,1,1)),:)=0;
        for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
           data.data(:,j)=zeros(1,length(data2.data(:,:)));
        end
        data2.data(:,:,i)=data.data(1:length(data2.data(:,1,1)),:);
    else
        for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
           data.data(:,j)=zeros(1,length(data2.data(:,:)));
        end
        data2.data(:,:,i)=data.data(1:length(data2.data(:,1,1)),:);
    end
end
 
end
 
% Input parameters
 
N_F = logdata.data(ramp_loc1,100); N_k = logdata.data(ramp_loc1,103);
F_min = logdata.data(ramp_loc1,98); F_max = logdata.data(ramp_loc1,99);
k_min = logdata.data(ramp_loc1,101); k_max = logdata.data(ramp_loc1,102);
k_sf = logdata.data(ramp_loc1,108);
alpha_init = logdata.data(ramp_loc1,106);
beta = logdata.data(ramp_loc1,107);
 
G_rand_init = spacedata.data(2,:);
Xc_rand_init = spacedata.data(3,:);
F_rand = spacedata.data(4,:);
k_rand = spacedata.data(5,:);
 
 
for i = 1:a
    alpha(i) = logdata.data(i,alpha_loc);
end
ramp = logdata.data(ramp_loc1,34)*100;
Fs = logdata.data(ramp_loc1,12);
pre = logdata.data(ramp_loc1,22);
pulse = logdata.data(ramp_loc1,23);
 
for i = 1:a
    for j = 118:(117+logdata.data(i,iter_loc))
        Xc_rand(i,j) = logdata.data(i,j);
    end
end
 
for i = 1:a
for j = (118+logdata.data(i,iter_loc)):(iter_loc+logdata.data(i,iter_loc)+logdata.data(i,iter_loc))
G_rand(i,j) = logdata.data(i,j);
end
end
 
% Iterate
 
for i =1:a
data0(:,:,i)=data2.data(:,:,i);
end
 
time = data2(:,1);
 
pre = pre*1e-3*Fs;
pulse = pulse*1e-3*Fs;
ramp = ramp*1e-3*Fs;
%pre = pre*1e-3*Fs;
%pulse = pulse*1e-3*Fs;
%ramp = ramp*1e-3*Fs;
 
for j = 1:a
    for i = 1:(logdata.data(j,iter_loc))
        Xc(:,i,j) = data0(:,(1+i),j);
        Xd(:,i,j) = data0(:,(1+logdata.data(j,iter_loc)+i),j);
        Delta(:,i,j) = data0(:,(1+2*logdata.data(j,iter_loc)+i),j);
        Gain(:,i,j) = data0(:,(1+3*logdata.data(j,iter_loc)+i),j);
    end
end
 
 %{
% Extract baselines and pulses
 
 
for j = 1:a
    for i = 1:(logdata.data(j,iter_loc))
        Xc_diff(:,i,j) = mean(Xc((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xc((pre/4):(pre/4*3),i,j));
        Xd_diff(:,i,j) = mean(Xd((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xd((pre/4):(pre/4*3),i,j));
        Delta_diff(:,i,j) = mean(Delta((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Delta((pre/4):(pre/4*3),i,j));
        Gain_diff(:,i,j) = mean(Gain((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Gain((pre/4):(pre/4*3),i,j));
    end
end  
 %}
offset =100*1e-3*Fs;
 
a_delay = logdata.data(ramp_loc1,15);
a_length = logdata.data(ramp_loc1,16);
a_freq = logdata.data(ramp_loc1,18);
a_amp = logdata.data(ramp_loc1,38);
 
 
for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
Xd_pulse(:,i,j) = Xd((pre+ceil(ramp)+offset):(pre+pulse-ceil(ramp)),i,j);
%Xd_pre(:,i,j) = Xd(((a_delay+a_length*1e-3*Fs+10):pre),i,j);
%Xd_post(:,i,j) = Xd((pre+ceil(ramp)+pulse+ceil(ramp)):length(Xd(:,i,j)),i,j);
%Xd_alpha(:,i,j) = Xd(a_delay*Fs/1000:(a_delay+a_length)*Fs/1000,i,j);
end
end
 
 
 
%% In-line alpha calibration
 
 
a_delay = logdata.data(ramp_loc1,15);
a_length = logdata.data(ramp_loc1,16);
a_freq = logdata.data(ramp_loc1,18);
a_amp = logdata.data(ramp_loc1,38);
 
alpha_sine = a_amp*sin(1e-3*2*a_freq*pi*time.data(a_delay*1e-3*Fs:(a_delay*1e-3*Fs+a_length*1e-3*Fs),1,1));
[w2(:,2),w2(:,1)] = pwelch(alpha_sine,[],[],1e6,Fs);
acal_amp = sqrt(w2(find(w2(:,1)==a_freq),2)*w2(find(w2(:,1)==a_freq),1))/2;
 
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
    acal_corr(:,i,j) = xcorr(alpha_sine,Xd_alpha(:,i,j),'unbiased');
    [w(:,2),w(:,1)] = pwelch(acal_corr(:,i,j),[],[],1e6,Fs);
    psd_w1(:,i,j) = w(:,1);
    psd_w2(:,i,j) = w(:,2);
    Xd_alpha_amp(i,j) = sqrt(w(find(w(:,1)==a_freq),2))/2;
    ascale(i,j) = acal_amp/Xd_alpha_amp(i,j);
end
end
 
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
    Xd_pulse(:,i,j)=Xd_pulse(:,i,j)*ascale(i,j);
    Xd_pre(:,i,j)=Xd_pre(:,i,j)*ascale(i,j);
    Xd_post(:,i,j)=Xd_post(:,i,j)*ascale(i,j);
end
end
 
 
 
 
%% Smooth and Calculate Amplitude and Frequency for Each
tic;
 
N=3;
 
min_freq = 5; % Fs/(length(Xd_pulse(:,1,1))/N);
max_freq = 100; %150;
max_freq2 = 45; %because intern alpha calib at 50 Hz right now.
 
  
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
Xd_pulse_z(:,i,j)=smooth(Xd_pulse(:,i,j),length(Xd_pulse(:,i,j))/N);
%Xd_pre_z(:,i,j)=smooth(Xd_pre(:,i,j),length(Xd_pre(:,i,j))/N);
%Xd_post_z(:,i,j)=smooth(Xd_post(:,i,j),length(Xd_post(:,i,j))/N);
end
end
 
Xd_pulse_center = Xd_pulse - Xd_pulse_z;
%Xd_pre_center = Xd_pre - Xd_pre_z;
%Xd_post_center = Xd_post - Xd_post_z;
 
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
Xd_RMS_amp(i,j)=std(Xd_pulse_center(:,i,j));
end
end
 
 
 
T=1/Fs;
L_p=length(Xd_pulse);
%L_b=length(Xd_pre);
NFFT_p=2^nextpow2(L_p);
%NFFT_b=2^nextpow2(L_b);
 
% Calculate power spectra
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
%q=Power_multi_traces(Xd_pulse_center(:,i,j),10,1/Fs,1);       % MTM              
%q2=Power_multi_traces(Xd_pre(:,i,j),10,1/Fs,1);
%[q(:,2),q(:,1)]=pwelch(Xd_pulse_center(:,i,j),[],[],1e5,Fs);   % pwelch
%[q2(:,2),q2(:,1)]=pwelch(Xd_pre(:,i,j),[],[],1e5,Fs);
%[q(:,2),q(:,1)]=spect(Xd_pulse_center(:,i,j),3,kais_beta,'ps',Fs);     % ASD
%[q2(:,2),q2(:,1)]=spect(Xd_pre_center(:,i,j),3,kais_beta,'ps',Fs);

h = fft(Xd_pulse_center(:,i,j),NFFT_p)/L_p;         % fft
qfft(:,2) = h(1:NFFT_p/2+1);
qfft(:,1) = Fs/2*linspace(0,1,NFFT_p/2+1);
%h2 = fft(Xd_pre_center(:,i,j),NFFT_b)/L_b;
%q2fft(:,2) = h2(1:NFFT_b/2+1);
%q2fft(:,1) = Fs/2*linspace(0,1,NFFT_b/2+1);

psd_q1fft(:,i,j)=qfft(:,1);    %frequency
psd_q2fft(:,i,j)=qfft(:,2);    %psd
%psd2_q1fft(:,i,j)=q2fft(:,1);
%psd2_q2fft(:,i,j)=q2fft(:,2);

h=psd(spectrum.periodogram('blackman'),Xd_pulse_center(:,i,j),'Fs',Fs);
q(:,1)=h.Frequencies;
q(:,2)=h.Data;
%h2=psd(spectrum.periodogram('blackman'),Xd_pre_center(:,i,j),'Fs',Fs);
%q2(:,1)=h2.Frequencies;
%q2(:,2)=h2.Data;
 
psd_q1(:,i,j)=q(:,1);    %frequency
psd_q2(:,i,j)=q(:,2);    %psd
%psd2_q1(:,i,j)=q2(:,1);
%psd2_q2(:,i,j)=q2(:,2);
p=find(q(:,2)==max(q(intersect(find(q(:,1)<max_freq),find(q(:,1)>min_freq)),2)));
%p2=find(q2(:,2)==max(q2(intersect(find(q(:,1)<max_freq2),find(q(:,1)>min_freq)),2)));
%p3=max(find(q2(:,1)<q(p,1)));

    Xd_freq(i,j)=q(p,1);
    Xd_amp(i,j)=sqrt(q(p,2)*(q(p+1,1)-q(p,1))*2); 


for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))                        %  pwelch method
[Pxx,Fxx]=pwelch(Xd_pulse_center(:,i,j),[],[],[],Fs);
psdwelch_q1(:,i,j)=Fxx;
psdwelch_q2(:,i,j)=Pxx;
p4=find(Pxx==max(Pxx(intersect(find(Fxx<max_freq),find(Fxx>min_freq)))));
Xd_amp_welch(i,j)=sqrt(Pxx(p4)*(Fxx(p4+1)-Fxx(p4))*2);
Xd_freq_welch(i,j)=Fxx(p4);
end
end
          % for PSD
%Xd_amp_base(i,j) = sqrt(q2(p3,2)*q2(p3,1));
%Xd_amp(i,j)=q(p,2);                         % for ASD
%Xd_amp_base(i,j)=q(p,1);
%Xd_amp(i,j)=2*abs(q(p,2));         % for fft
%Xd_amp_base(i,j)=2*abs(q(p,2));
%Xd_freq_base(i,j)=q2(p2,1);
%Xd_amp_ratio(i,j)=sqrt(q(p,2)*q(p,1))/sqrt(q2(p3,2)*q2(p3,1));
%Xd_freq_ratio(i,j)=q(p,1)/q2(p2,2);
end
end
 
%{
for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
q=Power_multi_traces(Xd_pulse_center(:,i,j),5,1/Fs,1);
q2=Power_multi_traces(Xd_pre(:,i,j),5,1/Fs,1);
psd_q1(:,i,j)=q(:,1);
psd_q2(:,i,j)=q(:,2);
psd2_q1(:,i,j)=q2(:,1);
psd2_q2(:,i,j)=q2(:,2);
p=find(q(:,2)==max(q(intersect(find(q(:,1)>min_freq),find(q(:,1)<max_freq)),2)));
p2=find(q2(:,2)==max(q2(intersect(find(q2(:,1)>min_freq),find(q2(:,1)<max_freq)),2)));
Xd_freq(i,j)=q(p,1);
Xd_freq_base(i,j)=q2(p2,1);
Xd_amp(i,j)=sqrt(q(p,2)*q(p,1));
Xd_amp_base(i,j)=sqrt(q2(p2,2)*q2(p2,1));
Xd_amp_base2(i,j)=sqrt(q2(p,2)*q2(p,1));
end
end
%}
 
toc
 
%% Autocorrelations

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
    pulse_corr(:,i,j)=xcorr(Xd_pulse_center(:,i,j),'biased');
 %   base_corr(:,i,j)=xcorr(Xd_pre_center(:,i,j),'biased');
    pulse_corr_env(:,i,j)=hilbert2(pulse_corr(:,i,j),Fs);
  %  base_corr_env(:,i,j)=hilbert2(base_corr(:,i,j),Fs);
    
end
end



%% Perform checks for bimodality and oscillations from probability distributions and power spectra
 
min_freq = Fs/(length(Xd_pulse(:,1,1))/N);
max_freq = 100;
dev = 3;
N = 2;


AIC = zeros(1,N); AICfit = zeros(1,N);
%NlogL = zeros(1,N);
clear cell;
obj = cell(1,N); fit1=cell(1,N);gof=cell(1,N);param=cell(1,N);

 %{
for j = 1:a
    for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
        prob_dist(:,i,j)=hist(Xd_pulse_center(:,i,j),100);
        [dip(i,j),p_dip(i,j),xlow_dip(i,j),xup(i,j)]=HartigansDipSignifTest(prob_dist(:,i,j),1000);
        %r = fit(psd_q1(300:1000,i,j),sqrt(psd_q1(300:1000,i,j).*psd_q2(300:1000,i,j)),'0*x+a','StartPoint',1);
        r.a = mean(sqrt(psd_q1(300:1000,i,j).*psd_q2(300:1000,i,j)));
        %r2 = fit(psd2_q1(300:1000,i,j),sqrt(psd2_q1(300:1000,i,j).*psd2_q2(300:1000,i,j)),'0*x+a','StartPoint',1);  
        r2.a = mean(sqrt(psd2_q1(300:1000,i,j).*psd2_q2(300:1000,i,j)));
        z = std(sqrt(psd_q1(300:1000,i,j).*psd_q2(300:1000,i,j)));
        z2 = std(sqrt(psd2_q1(300:1000,i,j).*psd2_q2(300:1000,i,j)));
   %     if Xd_amp(i,j) > r.a+dev*z && Xd_amp_base(i,j) > r2.a+dev*z2 &&  Xd_freq_base(i,j) >= min_freq && Xd_freq_base(i,j) <= max_freq && Xd_freq(i,j) >= min_freq && Xd_freq(i,j) <= max_freq
    %        Osc(i,j) = 1;       % Osc =  1 baseline and pulse oscillating
     %   elseif Xd_amp(i,j) > r.a+dev*z && Xd_amp_base(i,j) < r2.a+dev*z2 && Xd_freq(i,j) >= min_freq && Xd_freq(i,j) <= max_freq
      %      Osc(i,j) = 0.67;    % Osc = 0.67 baseline not oscillatin, pulse oscillating
      %  elseif Xd_amp(i,j) < r.a+dev*z && Xd_amp_base(i,j) > r2.a+dev*z2 && Xd_freq_base(i,j) >= min_freq && Xd_freq_base(i,j) <= max_freq 
      %      Osc(i,j) = 0.33;    % Osc = 0.33 baseline oscillating, pulse not oscillating 
      %  else
      %      Osc(i,j) = 0;       % Osc = 0 baseline and pulse not oscillating
      %  end
      if Xd_RMS_amp(i,j) > 9
          Osc(i,j) = 1;
      else
          Osc(i,j) = 0;
      end
end
end
 %}

 
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))                                            %%%% CHECK THIS LINE %%%%
    
[a2,b]=hist(Xd_pulse_center(:,i,j),200);


%g = fit(b',a',bimod_gauss,'StartPoint',[-20 20 0.5 10 10]);
%obj = gmdistribution.fit([b a]',2);
%c = pdf(obj,[b a]');
options = statset('Display','final');
%obj = gmdistribution.fit([b a]',2);
warning off

% The entire script iterates a number of times to account for variability
% in AIC, which can result in variable numbers of paramters. 

for q = 1:1

% Akaike information criterion (AIC) measures the goodness of fit. It is
% defined as 2k - 2ln(L), or two-times the number of parameters minus the
% maximum log-likelihood of the fit. This script iterates through normal
% distributions with different numbers of parameters (k=1:n). After the
% 'for' loop, the minimum AIC is used to determine the modality of this
% distribution.

for w = 1:N
    obj{w} = gmdistribution.fit(b',w,'Start','randSample');
%    AIC(k)= obj{k}.AIC;
%    NlogL(k) = obj{k}.NlogL;
    [fit1{w} gof{w} param{w}] = fit(b',a2',sprintf('%s%s', 'gauss',num2str(w)));
    AICfit(w) = 6*w + gof{w}.rmse;
end

%[minAIC,numComponents1] = min(AIC);
[minAICfit,numComponents1] = min(AICfit);
%[maxNlogL1,numComponents1] = max(NlogL);
numComponents(q)=numComponents1;
%NlogL(i)=maxNlogL1;

end

% Output model parameters
%model = obj{round(numComponents)}

Osc(i,j) = round(mean(numComponents));

warning on

end

Osc(i,j) = round(mean(numComponents));

end


for j = 1:a
    for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
        movingstd_psd(:,i,j)=sqrt(movingvar(psd2_q2(:,i,j),151));
        movingavg_psd(:,i,j)=moving(psd2_q2(:,i,j),151);
        
        movingstd_psd2(:,i,j)=sqrt(movingvar(psd2_q2(:,i,j),151));
        movingavg_psd2(:,i,j)=moving(psd2_q2(:,i,j),151);
%{
        for r=1:length(psd_q2(:,i,j))
            if psd_q2(r,i,j) > movingstd_psd(r,i,j)+movingavg_psd(r,i,j) 
                freq_peaks_ind(r)=psd_q1(r,i,j);
                amp_peaks_ind(r)=psd_q2(r,i,j);
            else
                freq_peaks_ind(r)=0;
                amp_peaks_ind(r)=0;                
            end 
    end
%}
        N_std = 1.5;
        
        freq_peaks_ind=psd_q1(find(psd_q2(1:500,i,j)>(movingstd_psd(1:500,i,j).*N_std+movingavg_psd(1:500,i,j))),i,j);     % 15:500 should sample first 100-150 Hz.... CHECK THIS
        amp_peaks_ind=psd_q2(find(psd_q2(1:500,i,j)>(movingstd_psd(1:500,i,j).*N_std+movingavg_psd(1:500,i,j))),i,j);
        if isempty(amp_peaks_ind)==1
            amp_peaks_ind=0;
            freq_peaks_ind=0;
        end
        amp_peaks(i,j)=sqrt(max(amp_peaks_ind)*freq_peaks_ind(find(max(amp_peaks_ind))));
        freq_peaks(i,j)=freq_peaks_ind(find(max(amp_peaks_ind)));
        clear freq_peaks_ind amp_peaks_ind
        
        freq_peaks_ind2=psd2_q1(find(psd2_q2(1:500,i,j)>(movingstd_psd2(1:500,i,j).*N_std+movingavg_psd2(1:500,i,j))),i,j);     % 15:500 should sample first 100-150 Hz.... CHECK THIS
        amp_peaks_ind2=psd2_q2(find(psd2_q2(1:500,i,j)>(movingstd_psd2(1:500,i,j).*N_std+movingavg_psd2(1:500,i,j))),i,j);
        amp_peaks2(i,j)=sqrt(max(amp_peaks_ind2)*freq_peaks_ind2(find(max(amp_peaks_ind2))));
        freq_peaks2(i,j)=freq_peaks_ind2(find(max(amp_peaks_ind2)));
        clear freq_peaks_ind2 amp_peaks_ind2        
    end
    
end


for j = 1:a
    for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
        if round(Xd_freq(i,j)) == round(freq_peaks(i,j)) | round(Xd_freq(i,j)) == round(freq_peaks2(i,j));
            Osc_p(i,j) = 1;
        else
            Osc_p(i,j) = 0;
    end
end
end

freq_peaks(amp_peaks<5)=0;

%% Eliminating noise
 
%Frequency
 
f_noise_max=100;
f_noise_min=Fs/(length(Xd_pulse(:,1,1))/N); %determined by the smoothing average : here, f_noise_min=1.6 Hz
 
%Xd_RMS_logdata.data(i+ramp_loc1-1(Xd_freq>f_noise_max)=0;
%Xd_RMS_amp(Xd_freq<f_noise_min)=0;
Xd_amp(Xd_freq>f_noise_max)=0;
Xd_amp(Xd_freq<f_noise_min)=0;
Xd_freq(Xd_freq>f_noise_max)=0;
Xd_freq(Xd_freq<f_noise_min)=0;
 
 
%Amplitude
 
%Xd_RMS_amp(Osc<0.5)=0;
Xd_amp(Osc<0.5)=0;
Xd_freq(Osc<0.5)=0;
 
Xd_freq(Xd_amp>100)=0;
Xd_amp(Xd_amp>100)=0;
 
 
 
 
 

%% Reshaping the vectors
 
 
%{ 
for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
Xd_amp_ratio(i,j)=Xd_amp(i,j)/Xd_amp_base(i,j); %%Xd_amp_base2 changed into 
Xd_freq_ratio(i,j)=Xd_freq(i,j)/Xd_freq_base(i,j);
end
end
 
 %}
 
 
 
k = size(Xd_RMS_amp);
 
Xd_RMS_amp = reshape(Xd_RMS_amp,1,k(1)*k(2));
Xd_freq = reshape(Xd_freq,1,k(1)*k(2));
%Xd_freq_ratio = reshape(Xd_freq_ratio,1,k(1)*k(2));
%Xd_amp_ratio = reshape(Xd_amp_ratio,1,k(1)*k(2));
Xd_amp = reshape(Xd_amp,1,k(1)*k(2));
%Gain_diff = reshape(Gain_diff,1,k(1)*k(2));
%Xc_diff = reshape(Xc_diff,1,k(1)*k(2));
%Osc = reshape(Osc,1,k(1)*k(2));
%freq_peaks = reshape(freq_peaks,1,k(1)*k(2));
%amp_peaks = reshape(amp_peaks,1,k(1)*k(2));
%freq_peaks2 = reshape(freq_peaks2,1,k(1)*k(2));
%amp_peaks2 = reshape(amp_peaks2,1,k(1)*k(2));
%Osc_p = reshape(Osc_p,1,k(1)*k(2));
Xd_amp_welch=reshape(Xd_amp_welch,1,k(1)*k(2));
Xd_freq_welch=reshape(Xd_freq_welch,1,k(1)*k(2));


%{
pulse_corr = reshape(pulse_corr,length(pulse_corr),1,k(1)*k(2));
pulse_corr_env = reshape(pulse_corr_env,length(pulse_corr_env),1,k(1)*k(2));
base_corr = reshape(base_corr,length(base_corr),1,k(1)*k(2));
base_corr_env = reshape(base_corr_env,length(base_corr_env),1,k(1)*k(2));
%} 
 
if length(k_rand)<(k(1)*k(2))
    for j = length(k_rand)+1:(k(1)*k(2))
        k_rand(j)=0;
    end
end
k_rand = k_rand(1:(k(1)*k(2)));
 
if length(F_rand)<(k(1)*k(2))
    for j = length(F_rand)+1:(k(1)*k(2))
        F_rand(j)=0;
    end
end
F_rand = F_rand(1:(k(1)*k(2)));
 
 
%%Ratio of RMS Amplitudes
 
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
Xd_RMS_amp_pre(i,j)=std(Xd_pre_center(:,i,j));
Xd_RMS_amp_post(i,j)=std(Xd_post_center(:,i,j));
Xd_RMS_amp_base(i,j)=mean([Xd_RMS_amp_pre(i,j) Xd_RMS_amp_post(i,j)]);
end
end
 
k = size(Xd_RMS_amp);
Xd_RMS_amp_base = reshape(Xd_RMS_amp_base,1,k(1)*k(2));
 
for j = 1:length(Xd_RMS_amp)
Xd_RMS_ratio(j)=Xd_RMS_amp(j)/Xd_RMS_amp_base(j);
end
 

 
 
%% Eliminating the zero points added for the last file that may not have as
%many points as the previous ones.
 
%%%%Displace this somewhere else%%%%%%
 
F_rand(find(Xd_RMS_amp==0))=[];
k_rand(find(Xd_RMS_amp==0))=[];
%Xc_diff(find(Xd_RMS_amp==0))=[];
%Gain_diff(find(Xd_RMS_amp==0))=[];
Xd_amp(find(Xd_RMS_amp==0))=[];
%Xd_amp_ratio(find(Xd_RMS_amp==0))=[];
%Xd_freq_ratio(find(Xd_RMS_amp==0))=[];
Xd_freq(find(Xd_RMS_amp==0))=[];
%Osc(find(Xd_RMS_amp==0))=[];
 
Xd_RMS_amp(find(Xd_RMS_amp==0))=[];
 
 
%% Plot with Interpolants

warning off;
 
figure(1);
F=TriScatteredInterp(k_rand',F_rand',Xd_RMS_amp');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.','MarkerSize',15);
colormap hot; colorbar;
title('RMS Amplitude (nm)'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
 
figure(2);
F=TriScatteredInterp(k_rand',F_rand',Xd_freq');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Major Frequency (Hz)'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
 
%{
figure(3);
F=TriScatteredInterp(k_rand',F_rand',Xc_diff');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none') hold on; axis tight;
%plot3(k_rand,F_rand,Xc_diff,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Command Displacement (nm)'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
 
figure(4);
F=TriScatteredInterp(k_rand',F_rand',Gain_diff');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none') hold on; axis tight;
%plot3(k_rand,F_rand,Gain_diff,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Gain'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
%}
 
figure(5);
F=TriScatteredInterp(k_rand',F_rand',Xd_amp');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%%plot3(k_rand,F_rand,Xd_amp,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Amplitude'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
 
figure(6);
F=TriScatteredInterp(k_rand',F_rand',Xd_amp_ratio');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_amp_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Amplitude Ratio'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
 
figure(7);
F=TriScatteredInterp(k_rand',F_rand',Xd_freq_ratio');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Frequency Ratio'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(8);
F=TriScatteredInterp(k_rand',F_rand',Osc');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Peaks in Histograms'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;


figure(9);
F=TriScatteredInterp(k_rand',F_rand',Osc_p');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Peaks from PSD (1=yes; 0=no)'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(10);
F=TriScatteredInterp(k_rand',F_rand',freq_peaks');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Frequency Peaks Blackman'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(11);
F=TriScatteredInterp(k_rand',F_rand',amp_peaks');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Amplitude Peaks Blackman'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;


figure(12);
F=TriScatteredInterp(k_rand',F_rand',amp_peaks2');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Amplitude Peaks FFT'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(13);
F=TriScatteredInterp(k_rand',F_rand',freq_peaks2');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz,'EdgeColor','none'); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_freq_ratio,'b.','MarkerSize',15);
colormap hot; colorbar;
title('Frequency Peaks FFT'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;



warning on;

 
%% Reordering as a matrix for imagesc
 
 
%C=zeros(N_F,N_k);
 
Ord_F=unique(F_rand);
Ord_k=unique(k_rand);
 
 
 
for i=1:N_F
        %p=find(F_rand==Ord_F(N_F+1-i));
        p=find(F_rand==Ord_F(i));
        for j=1:N_k-1
            q=find(k_rand(p)==Ord_k(j));
            Xd_RMS_amp_m(i,j)=Xd_RMS_amp(p(q));
            Xd_freq_m(i,j)=Xd_freq(p(q));
 %           Xc_diff_m(i,j)=Xc_diff(p(q));
 %           Gain_diff_m(i,j)=Gain_diff(p(q));
            Xd_amp_m(i,j)=Xd_amp(p(q));
 %           Xd_amp_ratio_m(i,j)=Xd_amp_ratio(p(q));
 %           Xd_freq_ratio_m(i,j)=Xd_freq_ratio(p(q));
 %           Osc_m(i,j)=Osc(p(q));
 %           Xd_amp_base_m(i,j)=Xd_amp_base(p(q));
 %           Xd_freq_base_m(i,j)=Xd_freq_base(p(q));
 %           Xd_RMS_amp_base_m(i,j)=Xd_RMS_amp_base(p(q));
            Xd_amp_welch_m(i,j)=Xd_amp_welch(p(q));
            Xd_freq_welch_m(i,j)=Xd_freq_welch(p(q));
        end
end
 
%Eliminate traces where baseline to different : 
%Figure out the mean+std of the base in amp and in frequency, and if the
%base is too deviant (ex above mean+3*std), not consider that point.
 
k1 = size(Xd_RMS_amp_m);
 %{
p=find(Xd_freq_base_m>(mean(reshape(Xd_freq_base_m,1,k1(1)*k1(2)))+3*std(reshape(Xd_freq_base_m,1,k1(1)*k1(2))))|Xd_freq_base_m<(mean(reshape(Xd_freq_base_m,1,k1(1)*k1(2)))-3*std(reshape(Xd_freq_base_m,1,k1(1)*k1(2)))));
 
Xd_freq_m(p)=NaN;
Xd_amp_m(p)=NaN;
Xd_amp_ratio_m(p)=NaN;
Xd_freq_ratio_m(p)=NaN;
 %}

load customcolormaps2
%Trace figures
 
%ty_m=sort(ty,'descend');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
 
figure(50); 
imagesc(tx*1e6,ty*1e12,Xd_RMS_amp_m)
colormap(red1); colorbar;
title('RMS Amplitude (nm)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(51);
imagesc(tx*1e6,ty*1e12,Xd_freq_m)
colormap(blue1); colorbar;
title('Major Frequency (Hz)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')

 
figure(54);
imagesc(tx*1e6,ty*1e12,Xd_amp_m);
colormap(red1); colorbar;
title('Amplitude'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')

figure(90); 
imagesc(tx*1e6,ty*1e12,Xd_RMS_amp_m)
colormap(gray2); colorbar;
title('RMS Amplitude (nm)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(91);
imagesc(tx*1e6,ty*1e12,Xd_freq_m)
colormap(gray2); colorbar;
title('Major Frequency (Hz)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')

 
figure(94);
imagesc(tx*1e6,ty*1e12,Xd_amp_m);
colormap(gray2); colorbar;
title('Amplitude'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
 
figure(52);
imagesc(tx*1e6,ty*1e12,Xd_amp_base_m);
colormap hot; colorbar;
title('Amplitude base'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(53);
imagesc(tx*1e6,ty*1e12,Xd_freq_base_m);
colormap hot; colorbar;
title('Frequency base'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')

figure(54);
imagesc(tx*1e6,ty*1e12,Xd_freq_welch_m);
colormap hot; colorbar;
title('Frequency (welch)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')

figure(55);
imagesc(tx*1e6,ty*1e12,Xd_amp_welch_m);
colormap hot; colorbar;
title('Amplitude (welch)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')



 
 
 
%{
figure(52);
imagesc(tx*1e6,ty*1e12,Xc_diff_m);
colormap hot; colorbar;
title('Command Displacement (nm)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(53);
imagesc(tx*1e6,ty*1e12,Gain_diff_m);
colormap gray; colorbar;
title('Gain'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
%}
 
figure(54);
imagesc(tx*1e6,ty*1e12,Xd_amp_m);
colormap hot; colorbar;
title('Amplitude'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(55);
imagesc(tx*1e6,ty*1e12,Xd_amp_ratio_m);
colormap hot; colorbar;
title('Amplitude Ratio'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(56);
imagesc(tx*1e6,ty*1e12,Xd_freq_ratio_m);
colormap hot; colorbar;
title('Frequency Ratio'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(57);
imagesc(tx*1e6,ty*1e12,Osc_m);
colormap gray; colorbar;
title('Oscillation map'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(58);
imagesc(tx*1e6,ty*1e12,Xd_RMS_amp_base_m);
colormap hot; colorbar;
title('Amplitude RMS base'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')

 
 
 
%% Taking the data and reordering it to get i=F and j=k
 
%F_rand and k_rand already in 1D vector
 
u=size(Xd_pulse_center);
 
Ord_F=unique(F_rand);
Ord_k=unique(k_rand);
 
 
 figure(60);
 
 
s=[logdata.data(1,iter_loc),a];
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(Xd_pulse_center(:,I(p),J(p)));
             axis([6e4 7e4 -40 40]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(Xd_pulse_center(:,I(p),J(p)));
              axis([0 1e4 -40 40]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
end
 
 
 
figure(61); 
 
 
s=[logdata.data(1,iter_loc),a];
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(Xd_pre_center(:,I(p),J(p)));
              axis([0 2e4 -60 60]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Baseline'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(Xd_pre_center(:,I(p),J(p)));
              axis([0 2e4 -40 40]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Baseline'])
       
      
    end
end
 
end
 
 
figure(62); 
 
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       hist(Xd_pulse_center(:,I(p),J(p)),100);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       hist(Xd_pulse_center(:,I(p),J(p)),100);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
end
 
figure(63); 
 
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(psd_q1(:,I(p),J(p)),psd_q2(:,I(p),J(p)));
             axis([0 30 0 50])
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(psd_q1(:,I(p),J(p)),psd_q2(:,I(p),J(p)));
             axis([0 30 0 50])
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
end

warning off;

figure(64); 
 
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(pulse_corr(length(pulse_corr)/2:length(pulse_corr),I(p),J(p)));
             axis([0 length(pulse_corr)/2 -100 100])
             hold on;
       plot(pulse_corr_env(length(pulse_corr_env)/2:length(pulse_corr_env),I(p),J(p)),'r');
             
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Pulse Correlation'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(pulse_corr(length(pulse_corr)/2:length(pulse_corr),I(p),J(p)));
             axis([0 length(base_corr)/2 -100 100])
             hold on;
       plot(pulse_corr_env(length(pulse_corr_env)/2:length(pulse_corr_env),I(p),J(p)),'r');             
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Pulse Correlation'])
       
      
    end
end
 
end
 



figure(65); 
 
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(base_corr(length(base_corr)/2:length(base_corr),I(p),J(p)));
             axis([0 length(base_corr)/2 -100 100])
             hold on;
       plot(base_corr_env(length(base_corr_env)/2:length(base_corr_env),I(p),J(p)),'r');
             
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Base Correlation'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(base_corr(length(base_corr)/2:length(base_corr),I(p),J(p)));
             axis([0 length(base_corr)/2 -100 100])
             hold on;
       plot(base_corr_env(length(base_corr_env)/2:length(base_corr_env),I(p),J(p)),'r');             
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Base Correlation'])
       
      
    end
end
 
end

warning on;

%Ideas : represent also the amplitude and frequency map of the pre-pulse
%time series. 


%%  Moving Variance
 
window = 5000;

for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
    movvar_Xd(:,i,j) = movingvar(Xd_pulse_center(:,i,j),window);
    movvar_post(:,i,j) = movingvar(Xd_post_center(:,i,j),window);
end
end


%% Modality of histograms


[a,b]=hist(Xd_pulse_center(:,1,2),200);


%g = fit(b',a',bimod_gauss,'StartPoint',[-20 20 0.5 10 10]);
%obj = gmdistribution.fit([b a]',2);
%c = pdf(obj,[b a]');
options = statset('Display','final');
%obj = gmdistribution.fit([b a]',2);
warning off

% The entire script iterates a number of times to account for variability
% in AIC, which can result in variable numbers of paramters. 

for i = 1:1

% Akaike information criterion (AIC) measures the goodness of fit. It is
% defined as 2k - 2ln(L), or two-times the number of parameters minus the
% maximum log-likelihood of the fit. This script iterates through normal
% distributions with different numbers of parameters (k=1:n). After the
% 'for' loop, the minimum AIC is used to determine the modality of this
% distribution.
N=2;
AIC = zeros(1,N);
%NlogL = zeros(1,N);
clear cell;
obj = cell(1,N); fit1=cell(1,N);gof=cell(1,N);param=cell(1,N);
for k = 1:N
    obj{k} = gmdistribution.fit(b',k,'Start','randSample');
%    AIC(k)= obj{k}.AIC;
%    NlogL(k) = obj{k}.NlogL;
    [fit1{k} gof{k} param{k}] = fit(b',a',sprintf('%s%s', 'gauss',num2str(k)));
    AICfit(k) = N*k - N*log(gof{k}.sse);
end

%[minAIC,numComponents1] = min(AIC);
[minAICfit,numcomponents1] = min(AICfit);
%[maxNlogL1,numComponents1] = max(NlogL);
numComponents(i)=numComponents1;
%NlogL(i)=maxNlogL1;

end

% Output model parameters
%model = obj{round(numComponents)}



warning on

%% Make movies
for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
    for k = 1:100
        plot(Xd_pulse_center(1:k*150,i,j));
        axis([0 7e4 -70 70])
        M(k,i,j)=getframe;
    end
end
end

for i = 1:7
for j = 1:7
movie2avi(M(:,i,j),sprintf('%s%s%s',num2str(i),num2str(j),'.avi'));
end
end
