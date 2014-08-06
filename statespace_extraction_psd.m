%% Import and Extract Data
% 
%   for x = 1:10
%       disp(x)
%   end
% 

clear all; close all;

% Initialization
iter_loc = 110;
alpha_loc=69;
ramp_loc=1;
path = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/';
date = '2012-11-18';
session = '.01';
ear = '1';
cell = '2';
user = 'JS';
sessions = ['123029 2';
'123139 3';
'123246 4';
'123327 5'];
statespace = '2';

spacefile = sprintf('%s%s%s%s %s%s %s%s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/State space',statespace,'.txt');
spacedata=importdata(spacefile,'\t',3);
logfile = sprintf('%s%s%s%s %s%s %s%s%s %s%s %s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',user,date,session,ear,cell,'.log');
logdata = importdata(logfile);

a = length(sessions(:,1));

for i = 1:a
file(i,:) = sprintf('%s%s%s%s %s%s %s%s%s %s%s %s %s %s %s%s',path,date,session,'/Ear',ear,'/Cell',cell,'/',user,date,session,ear,cell,'MS',sessions(i,:),'.txt');
end

for i = 1:a
data=importdata(file(i,:));

if i<a
    data2.data(:,:,i)=data.data;
else
    for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
        data.data(:,j)=zeros(1,length(data.data(:,:)));
    end
    data2.data(:,:,i)=data.data;
end

end

% Input parameters

N_F = logdata.data(ramp_loc,100); N_k = logdata.data(4,103);
F_min = logdata.data(ramp_loc,98); F_max = logdata.data(ramp_loc,99);
k_min = logdata.data(ramp_loc,101); k_max = logdata.data(ramp_loc,102);
k_sf = logdata.data(ramp_loc,108);
alpha_init = logdata.data(ramp_loc,106);
beta = logdata.data(ramp_loc,108);
%{
G_rand_init = spacedata.data(2,:);
Xc_rand_init = spacedata.data(3,:);
F_rand = spacedata.data(4,:);
k_rand = spacedata.data(5,:);
%}
G_rand_init = spacedata.data(:,2);
Xc_rand_init = spacedata.data(:,3);
F_rand = spacedata.data(:,4);
k_rand = spacedata.data(:,5);
F_rand=F_rand';
k_rand=k_rand';


for i = 1:a
    alpha(i) = logdata.data(i+ramp_loc-1,alpha_loc);
end
ramp = logdata.data(ramp_loc,34);
Fs = logdata.data(ramp_loc,11);
pre = logdata.data(ramp_loc,22);
pulse = logdata.data(ramp_loc,23);
a_delay = logdata.data(ramp_loc,15);
a_length = logdata.data(ramp_loc,16);

for i = 1:a
    for j = 121:(120+logdata.data(i+ramp_loc-1,iter_loc))
        Xc_rand(i,j) = logdata.data(i+ramp_loc-1,j);
    end
end

for i = 1:a
for j = (111+logdata.data(i+ramp_loc-1,iter_loc)):(iter_loc+logdata.data(i+ramp_loc-1,iter_loc)+logdata.data(i+ramp_loc-1,iter_loc))
G_rand(i,j) = logdata.data(i+ramp_loc-1,j);
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
a_delay = a_delay*1e-3*Fs;
a_length = a_length*1e-3*Fs;


for j = 1:a
    for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
        Xc(:,i,j) = data0(:,(1+i),j);
        Xd(:,i,j) = data0(:,(1+logdata.data(j+ramp_loc-1,iter_loc)+i),j);
        Delta(:,i,j) = data0(:,(1+2*logdata.data(j+ramp_loc-1,iter_loc)+i),j);
        Gain(:,i,j) = data0(:,(1+3*logdata.data(j+ramp_loc-1,iter_loc)+i),j);
    end
end


% Extract baselines and pulses


for j = 1:a
    for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
        Xc_diff(:,i,j) = mean(Xc((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xc((pre/4):(pre/4*3),i,j));
        Xd_diff(:,i,j) = mean(Xd((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xd((pre/4):(pre/4*3),i,j));
        Delta_diff(:,i,j) = mean(Delta((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Delta((pre/4):(pre/4*3),i,j));
        Gain_diff(:,i,j) = mean(Gain((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Gain((pre/4):(pre/4*3),i,j));
    end
end  

offset =0*1e-3*Fs;

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
Xd_pulse(:,i,j) = Xd((pre+ceil(ramp)+offset):(pre+pulse-ceil(ramp)),i,j);
Xd_pre(:,i,j) = Xd(((a_delay+a_length+10):pre),i,j);
Xd_post(:,i,j) = Xd((pre+ceil(ramp)+pulse+ceil(ramp)):length(Xd(:,i,j)),i,j);
Xd_alpha(:,i,j) = Xd(a_delay:(a_delay+a_length),i,j);
end
end

%% Smooth and Calculate Amplitude and Frequency for Each
tic;

N=3;

f_min=3;
f_max=200;

% Center traces about zero
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
Xd_pulse_z(:,i,j)=smooth(Xd_pulse(:,i,j),length(Xd_pulse(:,i,j))/N);
Xd_pre_z(:,i,j)=smooth(Xd_pre(:,i,j),length(Xd_pre(:,i,j))/N);
Xd_post_z(:,i,j)=smooth(Xd_post(:,i,j),length(Xd_post(:,i,j))/N);
end
end

Xd_pulse_center = Xd_pulse - Xd_pulse_z;
Xd_pre_center = Xd_pre - Xd_pre_z;
Xd_post_center = Xd_post - Xd_post_z;

%{
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
Xd_RMS_amp(i,j)=std(Xd_pulse_center(:,i,j));
end
end
%}

T=1/Fs;
L_p=length(Xd_pulse);
L_b=length(Xd_pre);
NFFT_p=2^nextpow2(L_p);
NFFT_b=2^nextpow2(L_b);

% Calculate power spectra
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
%q=Power_multi_traces(Xd_pulse_center(:,i,j),10,1/Fs,1);       % MTM              
%q2=Power_multi_traces(Xd_pre(:,i,j),10,1/Fs,1);
%[q(:,2),q(:,1)]=pwelch(Xd_pulse_center(:,i,j),[],[],1e5,Fs);   % pwelch
%[q2(:,2),q2(:,1)]=pwelch(Xd_pre(:,i,j),[],[],1e5,Fs);
%[q(:,2),q(:,1)]=spect(Xd_pulse_center(:,i,j),3,kais_beta,'ps',Fs);     % ASD
%[q2(:,2),q2(:,1)]=spect(Xd_pre_center(:,i,j),3,kais_beta,'ps',Fs);
%{
h = fft(Xd_pulse_center(:,i,j),NFFT_p)/L_p;         % fft
q(:,2) = h(1:NFFT_p/2+1);
q(:,1) = Fs/2*linspace(0,1,NFFT_p/2+1);
h2 = fft(Xd_pre_center(:,i,j),NFFT_b)/L_b;
q2(:,2) = h2(1:NFFT_p/2+1);
q2(:,1) = Fs/2*linspace(0,1,NFFT_b/2+1);
%}
h=psd(spectrum.periodogram('blackman'),Xd_pulse_center(:,i,j),'Fs',Fs);
q(:,1)=h.Frequencies;
q(:,2)=h.Data;
h2=psd(spectrum.periodogram('blackman'),Xd_pre_center(:,i,j),'Fs',Fs);
q2(:,1)=h2.Frequencies;
q2(:,2)=h2.Data;

psd_q1(:,i,j)=q(:,1);
psd_q2(:,i,j)=q(:,2);
psd2_q1(:,i,j)=q2(:,1);
psd2_q2(:,i,j)=q2(:,2);
p=find(q(:,2)==max(q(intersect(find(q(:,1)<f_max),find(q(:,1)>f_min)),2)));
p2=find(q2(:,2)==max(q2(intersect(find(q(:,1)<f_max),find(q(:,1)>f_min)),2)));
p3=find(max(q2(:,1)<q(p,1)));
Xd_freq(i,j)=q(p,1);
Xd_amp(i,j)=sqrt(q(p,2)*q(p,1));           % for PSD
Xd_amp_base(i,j) = sqrt(q2(p3,2)*q2(p3,1));
%Xd_amp(i,j)=q(p,2);                         % for ASD
%Xd_amp_base(i,j)=q(p,1);
%Xd_amp(i,j)=2*abs(q(p,2));         % for fft
%Xd_amp_base(i,j)=2*abs(q(p,2));
Xd_freq_base(i,j)=q2(p2,1);
Xd_amp_ratio(i,j)=sqrt(q(p,2)*q(p,1))/sqrt(q2(p3,2)*q2(p3,1));
Xd_freq_ratio(i,j)=q(p,1)/q2(p2,2);
end
end


toc

%% Perform checks for bimodality and oscillations from probability distributions and power spectra

min_freq = 3;
max_freq = 100;
dev = 3;

for j = 1:a
    for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
        prob_dist(:,i,j)=hist(Xd_pulse_center(:,i,j),100);
        [dip(i,j),p_dip(i,j),xlow_dip(i,j),xup(i,j)]=HartigansDipSignifTest(prob_dist(:,i,j),1000);
        r = fit(psd_q1(300:1000,i,j),sqrt(psd_q1(300:1000,i,j).*psd_q2(300:1000,i,j)),'0*x+a','StartPoint',1);
        r2 = fit(psd2_q1(300:1000,i,j),sqrt(psd2_q1(300:1000,i,j).*psd2_q2(300:1000,i,j)),'0*x+a','StartPoint',1);        
        z = sqrt(psd_q1(300:1000,i,j).*psd_q2(300:1000,i,j));
        if Xd_amp(i,j) >= r.a+dev*std(z) && Xd_amp_base(i,j) >= r2.a+dev*std(z)
            Osc(i,j) = 1;       % Osc =  1 (oscillating), 0.67 (oscillating, baseline not), 0.33 (baseline oscillating, pulse not) 0 (baseline also not oscillating)
        elseif Xd_amp(i,j) >= r.a+dev*std(z) && Xd_amp_base(i,j) < r2.a+dev*std(z)
            Osc(i,j) = 0.67;
        elseif Xd_amp_base(i,j) >= r2.a+std(z) && Xd_freq_base(i,j) >= min_freq && Xd_freq_base(i,j) <= max_freq && Xd_freq(i,j) >= min_freq && Xd_freq(i,j) <= max_freq
            Osc(i,j) = 0.33;
        else
            Osc(i,j) = 0;
        end
    end
end


        %% Plot Frequency and Amplitude Maps

k = size(Xd_freq);

%Xd_RMS_amp = reshape(Xd_RMS_amp,1,k(1)*k(2));

%Xd_psd_amp = reshape(Xd_psd_amp,1,k(1)*k(2));
Xd_freq = reshape(Xd_freq,1,k(1)*k(2));
Xd_freq_ratio = reshape(Xd_freq_ratio,1,k(1)*k(2));
Xd_amp_ratio = reshape(Xd_amp_ratio,1,k(1)*k(2));
Xd_amp = reshape(Xd_amp,1,k(1)*k(2));
Gain_diff = reshape(Gain_diff,1,k(1)*k(2));
Xc_diff = reshape(Xc_diff,1,k(1)*k(2));
dip = reshape(dip,1,k(1)*k(2));
p_dip = reshape(p_dip,1,k(1)*k(2));
Osc = reshape(Osc,1,k(1)*k(2));

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
tri = delaunay(k_rand,F_rand);
%{
figure(1); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_RMS_amp);
axis vis3d; colormap hot;
title('RMS Amplitude (nm)'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;
%}
figure(2); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_freq);
axis vis3d; colormap hot;
title('Major Frequency (Hz)'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;

figure(3); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xc_diff);
axis vis3d; colormap hot;
title('Command Displacement (nm)'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;

figure(4); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Gain_diff);
axis vis3d; colormap hot;
title('Gain'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;

figure(5); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_amp);
axis vis3d; colormap hot;
title('Amplitude'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;

figure(6); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_amp_ratio);
axis vis3d; colormap hot;
title('Power Ratio'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;

figure(7); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_freq_ratio);
axis vis3d; colormap hot;
title('Frequency Ratio'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;
%{
figure(8); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_psd_amp);
axis vis3d; colormap hot;
title('PSD Amplitude'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;
%}
figure(9); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,dip);
axis vis3d; colormap hot;
title('Dip Test for Bimodality'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;

figure(10); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,p_dip);
axis vis3d; colormap hot;
title('Dip p-value'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;
%% Ratio of RMS Amplitudes
%{
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
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
%}
%% Plot old and new RMS Amplitudes
%{
figure(1); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_RMS_amp);
axis vis3d; colormap hot;
title('RMS Amplitude (nm)'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;

figure(5); hold on;
h = trisurf(tri,k_rand*1e6,F_rand*1e12,Xd_RMS_ratio);
axis vis3d; colormap hot;
title('RMS Amplitude Ratio (pulse/baseline)'); xlabel('stiffness (uN/m)'); ylabel('force (pN)');
hold off;
%}
%% Plot with Interpolants
%{
figure(1);
F=TriScatteredInterp(k_rand',F_rand',Xd_RMS_amp');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('RMS Amplitude (nm)'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
%}
figure(2);
F=TriScatteredInterp(k_rand',F_rand',Xd_freq');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Major Frequency (Hz)'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(3);
F=TriScatteredInterp(k_rand',F_rand',Xc_diff');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Command Displacement (nm)'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(4);
F=TriScatteredInterp(k_rand',F_rand',Gain_diff');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Gain'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(5);
F=TriScatteredInterp(k_rand',F_rand',Xd_amp');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Amplitude'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(6);
F=TriScatteredInterp(k_rand',F_rand',Xd_amp_ratio');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Amplitude Ratio'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(7);
F=TriScatteredInterp(k_rand',F_rand',Xd_freq_ratio');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Frequency Ratio'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
%{
figure(8);
F=TriScatteredInterp(k_rand',F_rand',Xd_psd_amp');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('PSD Amplitude'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;
%}
figure(9);
F=TriScatteredInterp(k_rand',F_rand',dip');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Dip test for bimodality'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(10);
F=TriScatteredInterp(k_rand',F_rand',p_dip');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('p-value for dip test'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

figure(11);
F=TriScatteredInterp(k_rand',F_rand',Osc');
tx=min(k_rand):(max(k_rand)-min(k_rand))/N_k:max(k_rand);
ty=min(F_rand):(max(F_rand)-min(F_rand))/N_F:max(F_rand);
[qx,qy]=meshgrid(tx,ty);
qz=F(qx,qy);
surf(qx,qy,qz); hold on; axis tight;
%plot3(k_rand,F_rand,Xd_RMS_amp,'b.');
colormap hot; colorbar;
title('Oscillating'); xlabel('stiffness (N/m)'); ylabel('force (N)');
hold off;

%% Plot data for visualization
u=size(Xd_pulse_center);

dim1=3;
dim2=3;

for j=1:a
    for i = 1:logdata.data(j+ramp_loc-1,iter_loc)
       figure(10+j);
       subplot(dim1,dim2,i)
       hold on;
       hist(Xd_pulse_center(:,i,j),100);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],i,j)),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],i,j)),3), ' dip=',num2str(dip(sub2ind([u(2) u(3)],i,j)),3), 'p_dip=',num2str(p_dip(sub2ind([u(2) u(3)],i,j)),3)])
    end
end

for j=1:a
    for i = 1:logdata.data(j+ramp_loc-1,iter_loc)
       figure(10+j+a);
       subplot(dim1,dim2,i)
       hold on;
       plot(Xd_pulse_center(:,i,j));
              axis([0 2e4 -40 40]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],i,j)),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],i,j)),3), 'PULSE'])
    end
end


for j=1:a
    for i = 1:logdata.data(j+ramp_loc-1,iter_loc)
       figure(10+j+2*a);
       subplot(dim1,dim2,i)
       hold on;
       plot(Xd_pre_center(:,i,j));
       axis([0 2e4 -40 40]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],i,j)),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],i,j)),3), 'BASELINE'])
    end
end

for j=1:a
    for i = 1:logdata.data(j+ramp_loc-1,iter_loc)
       figure(10+j+3*a);
       subplot(dim1,dim2,i)
       hold on;
       plot(psd_q1(:,i,j),psd_q2(:,i,j));
       axis([0 15 0 50])
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],i,j)),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],i,j)),3), 'PSD'])
    end
end

%% Remove zero values

F_rand(find(Xd_freq==0))=[];
k_rand(find(Xd_freq==0))=[];
Xc_diff(find(Xd_freq==0))=[];
Gain_diff(find(Xd_freq==0))=[];
Xd_amp(find(Xd_freq==0))=[];
Xd_amp_ratio(find(Xd_freq==0))=[];
Xd_freq_ratio(find(Xd_freq==0))=[];
Xd_freq(find(Xd_freq==0))=[];
Xd_RMS_amp(find(Xd_freq==0))=[];
Osc(find(Xd_freq==0))=[];

%% Exclude extreme values

 
%Frequency
 
f_noise_max=100;
f_noise_min=Fs/(length(Xd_pulse(:,1,1))/N); %determined by the smoothing average : here, f_noise_min=1.6 Hz
r1=find(Xd_freq>f_noise_max); 
r2= find(Xd_freq<f_noise_min);
Xd_freq(r1)=0;
Xd_freq(Xd_freq<f_noise_min)=0;

Xd_freq_ratio(r1)=0;
Xd_freq_ratio(r2)=0;

%F_rand(r1)=0;
%F_rand(r2)=0;

%k_rand(r1)=0;
%k_rand(r2)=0;

Xd_amp(r1)=0;
Xd_amp(r2)=0;

Xd_amp_ratio(r1)=0;
Xd_amp_ratio(r2)=0;

Xc_diff(r1)=0;
Xc_diff(r2)=0;

Gain_diff(r1)=0;
Gain_diff(r2)=0;

%% Reordering as a matrix for imagesc
 
 
%C=zeros(N_F,N_k);
 
Ord_F=unique(F_rand);
Ord_k=unique(k_rand);

for i=1:N_F
        %p=find(F_rand==Ord_F(N_F+1-i));
        p=find(F_rand==Ord_F(i));
        for j=1:N_k
            q=find(k_rand(p)==Ord_k(j));
    %        Xd_RMS_amp_m(i,j)=Xd_RMS_amp(p(q));
            Xd_freq_m(i,j)=Xd_freq(p(q));
            Xc_diff_m(i,j)=Xc_diff(p(q));
            Gain_diff_m(i,j)=Gain_diff(p(q));
            Xd_amp_m(i,j)=Xd_amp(p(q));
            Xd_amp_ratio_m(i,j)=Xd_amp_ratio(p(q));
            Xd_freq_ratio_m(i,j)=Xd_freq_ratio(p(q));
            Osc_m(i,j)=Osc(p(q));

        end
end
 
%tx=Ord_k;
%ty=Ord_F;
%ty_m=sort(ty,'descend');
 %{
figure(8); 
imagesc(tx*1e6,ty*1e12,Xd_RMS_amp_m)
colormap hot; colorbar;
title('RMS Amplitude (nm)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 %}

figure(9);
imagesc(tx*1e6,ty*1e12,Xd_freq_m)
colormap bone; colorbar;
title('Major Frequency (Hz)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(10);
imagesc(tx*1e6,ty*1e12,Xc_diff_m);
colormap cool; colorbar;
title('Command Displacement (nm)'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(11);
imagesc(tx*1e6,ty*1e12,Gain_diff_m);
colormap cool; colorbar;
title('Gain'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(12);
imagesc(tx*1e6,ty*1e12,Xd_amp_m);
colormap bone; colorbar;
title('Amplitude'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(13);
imagesc(tx*1e6,ty*1e12,Xd_amp_ratio_m);
colormap hot; colorbar;
title('Amplitude Ratio'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 
figure(14);
imagesc(tx*1e6,ty*1e12,Xd_freq_ratio_m);
colormap hot; colorbar;
title('Frequency Ratio'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')

figure(15);
imagesc(tx*1e6,ty*1e12,Osc_m);
colormap jet; colorbar;
title('Oscillating'); xlabel('stiffness (µN/m)'); ylabel('force (pN)');
set(gca,'YDir','normal')
 

%% Taking the data and reordering it to get i=F and j=k
 
%F_rand and k_rand already in 1D vector
 
u=size(Xd_pulse_center);
 
Ord_F=unique(F_rand);
Ord_k=unique(k_rand);
 
 
 figure(40);
 
 
s=[logdata.data(1,iter_loc),a];
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(Xd_pulse_center(:,I(p),J(p)));
              axis([0 2e4 -40 40]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
 
 
%{
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(Xd_pulse_center(:,I(p),J(p)));
              axis([1e4 1.5e4 -20 20]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
    end
end
 
 
%}

%% In-line alpha calibration

a_delay = logdata.data(ramp_loc,15);
a_length = logdata.data(ramp_loc,16);
a_freq = logdata.data(ramp_loc,18);
a_amp = logdata.data(ramp_loc,38);

alpha_sine = a_amp*sin(1e-3*2*a_freq*pi*time.data(a_delay*1e-3*Fs:(a_delay*1e-3*Fs+a_length*1e-3*Fs),1,1));
[w2(:,2),w2(:,1)] = pwelch(alpha_sine,[],[],1e6,Fs);
acal_amp = sqrt(w2(find(w2(:,1)==a_freq),2)*w2(find(w2(:,1)==a_freq),1))/2;

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
    acal_corr(:,i,j) = xcorr(alpha_sine,Xd_alpha(:,i,j),'unbiased');
    [w(:,2),w(:,1)] = pwelch(acal_corr(:,i,j),[],[],1e6,Fs);
    psd_w1(:,i,j) = w(:,1);
    psd_w2(:,i,j) = w(:,2);
    Xd_alpha_amp(i,j) = sqrt(w(find(w(:,1)==a_freq),2))/2;
    ascale(i,j) = acal_amp/Xd_alpha_amp(i,j);
end
end

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc-1,iter_loc))
    Xd_pulse(:,i,j)=Xd_pulse(:,i,j)*ascale(i,j);
    Xd_pre(:,i,j)=Xd_pre(:,i,j)*ascale(i,j);
    Xd_post(:,i,j)=Xd_post(:,i,j)*ascale(i,j);
end
end


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
              axis([0 2e4 -40 40]);
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
              axis([0 2e4 -40 40]);
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
              axis([0 2e4 -40 40]);
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
       plot(Xd_pre_center(:,I(p),J(p)));
              axis([0 2e4 -40 40]);
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'PULSE'])
       
      
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


