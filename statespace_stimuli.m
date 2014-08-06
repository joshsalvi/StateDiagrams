%% Extract Data


clear all; close all;
 
iter_loc = 110; %logdata.data(i,iter_loc) : number of state points in the file i
alpha_loc=69;
ramp_loc1=1; %ramp_loc2=37;
path = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/';
date = '2014-03-19';
session = '.01';
ear = '1';
cell = '9';
user = 'JS';
sessions = ['135014 2';
'135139 3';
'135301 4'];
statespace = '1';
 
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
ramp = logdata.data(ramp_loc1,34);
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
    for i = 1:(logdata.data(ramp_loc1,iter_loc))
        Xc(:,i,j) = data0(:,(1+i),j);
        Xd(:,i,j) = data0(:,(1+logdata.data(j,iter_loc)+i),j);
       Delta(:,i,j) = data0(:,(1+2*logdata.data(j,iter_loc)+i),j);
%        Gain(:,i,j) = data0(:,(1+3*logdata.data(j,iter_loc)+i),j);
    end
end


%{
 
for j = 1:a
    for i = 1:(logdata.data(ramp_loc1,iter_loc))
        Xc_diff(:,i,j) = mean(Xc((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xc((pre/4):(pre/4*3),i,j));
        Xd_diff(:,i,j) = mean(Xd((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xd((pre/4):(pre/4*3),i,j));
   %     Delta_diff(:,i,j) = mean(Delta((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Delta((pre/4):(pre/4*3),i,j));
   %     Gain_diff(:,i,j) = mean(Gain((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Gain((pre/4):(pre/4*3),i,j));
    end
end  
%} 
offset =100*1e-3*Fs;
 
a_delay = logdata.data(ramp_loc1,15);
a_length = logdata.data(ramp_loc1,16);
a_freq = logdata.data(ramp_loc1,18);
a_amp = logdata.data(ramp_loc1,38);
 
 
for j = 1:a
for i = 1:(logdata.data(ramp_loc1,iter_loc))
Xd_pulse(:,i,j) = Xd((pre+ceil(ramp)+offset):(pre+pulse-ceil(ramp)),i,j);
Xd_pre(:,i,j) = Xd(((a_delay+a_length*1e-3*Fs+10):pre),i,j);
Xd_post(:,i,j) = Xd((pre+ceil(ramp)+pulse+ceil(ramp)):length(Xd(:,i,j)),i,j);
Xd_alpha(:,i,j) = Xd(a_delay*Fs/1000:(a_delay+a_length)*Fs/1000,i,j);
Xc_pulse(:,i,j) = Xc((pre+ceil(ramp)+offset):(pre+pulse-ceil(ramp)),i,j);
Delta_pulse(:,i,j) = Delta((pre+ceil(ramp)+offset):(pre+pulse-ceil(ramp)),i,j);
end
end



%% In-line alpha calibration
 
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
 
 

%% Extract stimulus delivery

choice = menu('Stimulus Type','DC (none)','Pulse','Discrete Frequencies','Frequency Sweep','Ramp');

switch choice
    case 1          % DC (none)    
        disp('Selected DC');
        
         for j = 1:a 
            for i = 1:(logdata.data(j,iter_loc))
                Xc_diff(:,i,j) = mean(Xc((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xc((pre/4):(pre/4*3),i,j));
                Xd_diff(:,i,j) = mean(Xd((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Xd((pre/4):(pre/4*3),i,j));
                Delta_diff(:,i,j) = mean(Delta((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Delta((pre/4):(pre/4*3),i,j));
                Gain_diff(:,i,j) = mean(Gain((pre+(pulse/4)):(pre+(pulse/2)*1.5),i,j))-mean(Gain((pre/4):(pre/4*3),i,j));
            end
        end 
    
    
    case 2          % Pulse
        disp('Selected PULSE');
        
        delay1=logdata.data(1,114)*Fs*1e-3;
        pulseduration=logdata.data(1,115)*Fs*1e-3;
        pulseamp=logdata.data(1,117)*Fs*1e-3;
        totalperiod=logdata.data(1,119)*Fs*1e-3;
        os = totalperiod-length(Xd_pulse);
        
        for j = 1:a
            for i = 1:(logdata.data(j,iter_loc))
                Xd_pulse_prepulse(:,i,j) = Xd_pulse(1:(delay1-os),i,j);
                Xd_pulse_pulse(:,i,j) = Xd_pulse((delay1-os):((delay1-os)+pulseduration),i,j);
                Xd_pulse_postpulse(:,i,j) = Xd_pulse((delay1-os)+pulseduration:(totalperiod-os),i,j);
            end
        end
        
        
    case 3          % Discrete Frequencies
        disp('Selected DISCRETE FREQUENCIES');
        
        delay1=logdata.data(1,114)*Fs*1e-3;

        
        choice2=menu('Sweep Type','Linear','Array','Exponential','Amplitude Array');
        
        switch choice2
            case 1
                disp('Linear');
                
                delay1=logdata.data(1,114)*Fs*1e-3;
                totalperiod = logdata.data(1,119)*Fs*1e-3;
                
                cycles=logdata.data(1,116);
                freq1=logdata.data(1,117);
                freq2=logdata.data(1,118);
                numfreq=logdata.data(1,123);
                freqinc = (freq2-freq1)/(numfreq-1);        % Linear

                offset = 750;           % If Xd_pulse is cut off from the beginning, set this value accordingly.
                
                for i = 1:numfreq
                    freq_stim(i)=freq1+freqinc*(i-1);       % Linear array of frequencies
                    stim_time(i)=cycles/freq_stim(i)*Fs;   % Duration defined by CYCLES
                end
                
                
                for j = 1:a
                for i = 1:(logdata.data(j,iter_loc))
                    for p = 1:numfreq
                    if p>1
                        clear stim_pad;
                        stim_pad = Xd_pulse((delay1*p+sum(stim_time(1:p-1)))-offset:(delay1*p+sum(stim_time(1:p)))-offset,i,j);
                        stim_pad(length(stim(:,1,1,1)))=0;
                        stim(:,p,i,j)=stim_pad;
                    else
                        stim(:,p,i,j)=Xd_pulse((delay1*p)-offset:(delay1*p+stim_time(p)-offset),i,j);
                    end
                end
                end
                end
                    
                
            case 2
                disp('Array');
                
                delay1=logdata.data(1,114)*Fs*1e-3;
                totalperiod = logdata.data(1,119)*Fs*1e-3;
                
                cycles=logdata.data(1,116);
                freq1=logdata.data(1,117);
                freq2=logdata.data(1,118);
                numfreq=logdata.data(1,123);
                
                freq_stim=[5	9	13	15	17	19	21	23	25	27	30	33	36	38	40	43	46	50	55	60	70	80];     %%%% INPUT FREQUENCY ARRAY %%%%
                numfreq=length(freq_stim);
                for i = 1:numfreq
                     stim_time(i)=cycles/freq_stim(i)*Fs;   % Duration defined by CYCLES
                end
                
                
                
                for j = 1:a
                for i = 1:(logdata.data(j,iter_loc))
                    for p = 1:numfreq
                    if p>1
                        clear stim_pad;
                        stim_pad = Xd_pulse((delay1*p+sum(stim_time(1:p-1)))-offset:(delay1*p+sum(stim_time(1:p)))-offset,i,j);
                        input_pad = Xc_pulse((delay1*p+sum(stim_time(1:p-1)))-offset:(delay1*p+sum(stim_time(1:p)))-offset,i,j);
                        Delta_pad = Delta_pulse((delay1*p+sum(stim_time(1:p-1)))-offset:(delay1*p+sum(stim_time(1:p)))-offset,i,j);
                        stim_pad(length(stim(:,1,1,1)))=0;
                        input_pad(length(stim(:,1,1,1)))=0;
                        Delta_pad(length(stim(:,1,1,1)))=0;
                        stim(:,p,i,j)=stim_pad;
                        input(:,p,i,j)=input_pad;
                        Delta_in(:,p,i,j)=Delta_pad;
                    else
                        stim(:,p,i,j)=Xd_pulse((delay1*p)-offset:(delay1*p+stim_time(p)-offset),i,j);
                        input(:,p,i,j)=Xc_pulse((delay1*p)-offset:(delay1*p+stim_time(p)-offset),i,j);
                        Delta_in(:,p,i,j)=Delta_pulse((delay1*p)-offset:(delay1*p+stim_time(p)-offset),i,j);
                    end
                end
                end
                end

            case 3
                disp('Exponential');
                
            case 4
                disp('Amplitude Sweep, ARRAY');
                
                delay1=logdata.data(1,114)*Fs*1e-3;
                totalperiod = logdata.data(1,119)*Fs*1e-3;
                
                cycles=logdata.data(1,116);
                freq1=logdata.data(1,117);
                freq2=logdata.data(1,118);
                %numfreq=logdata.data(1,123);
                
                freq_stim=freq1;                                             %%%% INPUT FREQUENCY %%%%
                amp_stim=[5	10	15	20	25	30	35	40	50	60	70	80	90	100	150	200];     %%%% INPUT AMPLITUDE ARRAY %%%%
                numfreq=length(amp_stim);
                
                for i = 1:numfreq
                     stim_time(i)=cycles/freq_stim*Fs;   % Duration defined by CYCLES
                end
                
                
                
                for j = 1:a
                for i = 1:(logdata.data(j,iter_loc))
                    for p = 1:numfreq
                    if p>1
                        clear stim_pad;
                        stim_pad = Xd_pulse((delay1*p+sum(stim_time(1:p-1)))-offset:(delay1*p+sum(stim_time(1:p)))-offset,i,j);
                        input_pad = Xc_pulse((delay1*p+sum(stim_time(1:p-1)))-offset:(delay1*p+sum(stim_time(1:p)))-offset,i,j);
                        Delta_pad = Delta_pulse((delay1*p+sum(stim_time(1:p-1)))-offset:(delay1*p+sum(stim_time(1:p)))-offset,i,j);
                        stim_pad(length(stim(:,1,1,1)))=0;
                        input_pad(length(stim(:,1,1,1)))=0;
                        Delta_pad(length(stim(:,1,1,1)))=0;
                        stim(:,p,i,j)=stim_pad;
                        input(:,p,i,j)=input_pad;
                        Delta_in(:,p,i,j)=Delta_pad;
                    else
                        stim(:,p,i,j)=Xd_pulse((delay1*p)-offset:(delay1*p+stim_time(p)-offset),i,j);
                        input(:,p,i,j)=Xc_pulse((delay1*p)-offset:(delay1*p+stim_time(p)-offset),i,j);
                        Delta_in(:,p,i,j)=Delta_pulse((delay1*p)-offset:(delay1*p+stim_time(p)-offset),i,j);
                    end
                end
                end
                end

                
                
        end
        
    case 4          % Frequency Sweep
        disp('Selected FREQUENCY SWEEP');
        
    case 5          % Ramp
        disp('Selected RAMP');
        
end
 
%% Case 1 - DC 

if choice == 1
    
    
    
    
    
    
    
    
    
else
    disp('Incorrect Choice');
end

%% Case 2 - Pulse

if choice == 2
    
N=1;
 
min_freq = 5; % Fs/(length(Xd_pulse(:,1,1))/N);
max_freq = 100; %150;
max_freq2 = 45; %because intern alpha calib at 50 Hz right now.
 
 
for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
    Xd_pulse_prepulse_z(:,i,j)=smooth(Xd_pulse_prepulse(:,i,j),length(Xd_pulse_prepulse(:,i,j))/N);
    Xd_pulse_pulse_z(:,i,j)=smooth(Xd_pulse_pulse(:,i,j),length(Xd_pulse_pulse(:,i,j))/N);
    Xd_pulse_postpulse_z(:,i,j)=smooth(Xd_pulse_postpulse(:,i,j),length(Xd_pulse_postpulse(:,i,j))/N);
end
end

Xd_pulse_prepulse_center = Xd_pulse_prepulse - Xd_pulse_prepulse_z;
Xd_pulse_pulse_center = Xd_pulse_pulse - Xd_pulse_pulse_z;
Xd_pulse_postpulse_center = Xd_pulse_postpulse - Xd_pulse_postpulse_z;

for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
Xd_pulse_prepulse_RMS_amp(i,j)=std(Xd_pulse_prepulse_center(round(length(Xd_pulse_prepulse)/8):round(length(Xd_pulse_prepulse)-length(Xd_pulse_prepulse)/8),i,j));
Xd_pulse_postpulse_RMS_amp(i,j)=std(Xd_pulse_postpulse_center(round(length(Xd_pulse_postpulse)/8):round(length(Xd_pulse_postpulse)-length(Xd_pulse_postpulse)/8),i,j));
end
end

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
h_pre=psd(spectrum.periodogram('blackman'),Xd_pulse_prepulse_center(:,i,j),'Fs',Fs);
q_pre(:,1)=h_pre.Frequencies;
q_pre(:,2)=h_pre.Data;
h_post=psd(spectrum.periodogram('blackman'),Xd_pulse_postpulse_center(:,i,j),'Fs',Fs);
q_post(:,1)=h_post.Frequencies;
q_post(:,2)=h_post.Data;

psd_q1_pre(:,i,j)=q_pre(:,1);    %frequency PRE
psd_q2_pre(:,i,j)=q_pre(:,2);    %psd PRE
psd_q1_post(:,i,j)=q_post(:,1);    %frequency POST
psd_q2_post(:,i,j)=q_post(:,2);    %psd POST

p_pre=find(q_pre(:,2)==max(q_pre(intersect(find(q_pre(:,1)<max_freq),find(q_pre(:,1)>min_freq)),2)));
p_post=find(q_post(:,2)==max(q_post(intersect(find(q_post(:,1)<max_freq),find(q_post(:,1)>min_freq)),2)));

Xd_freq_pre(i,j)=q_pre(p_pre,1);
Xd_amp_pre(i,j)=sqrt(q_pre(p_pre,2)*q_pre(p_pre,1)); 
Xd_freq_post(i,j)=q_post(p_post,1);
Xd_amp_post(i,j)=sqrt(q_post(p_post,2)*q_pre(p_post,1)); 

end
end

    
% Exponential fit

N=30;
start=300;      % If Xd_pulse is cut off at the beginnning, set this value accordingly

t = linspace(0,length(Xd_pulse_postpulse(start:round(length(Xd_pulse_postpulse)/N),1,1))/Fs,length(Xd_pulse_postpulse(start:round(length(Xd_pulse_postpulse)/N),1,1)));
warning on
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
    [post_fit gof]=fit(t',-1*Xd_pulse_postpulse(start:round(length(Xd_pulse_postpulse)/N),i,j),'exp2');
    if post_fit.b < 0  
        post_fit_a(i,j)=post_fit.a;
        post_fit_b(i,j)=post_fit.b;
    else if post_fit.d < 0
        post_fit_a(i,j)=post_fit.a;
        post_fit_b(i,j)=post_fit.b;
    else
        post_fit_a(i,j)=post_fit.a;
        post_fit_b(i,j)=post_fit.b;
        end
    end
    
    post_rsquare(i,j) = gof.rsquare;
    exp_post(:,i,j)=post_fit.a*exp(post_fit.b*t);

end
end
warning off    


% Reorder matrices

Ord_F=unique(F_rand);
Ord_k=unique(k_rand);

for i=1:N_F
        %p=find(F_rand==Ord_F(N_F+1-i));
        p=find(F_rand==Ord_F(i));
        for j=1:N_k
            q=find(k_rand(p)==Ord_k(j));
            Xd_freq_pre_m(i,j)=Xd_freq_pre(p(q));
            Xd_freq_post_m(i,j)=Xd_freq_post(p(q));
            Xd_amp_pre_m(i,j)=Xd_amp_pre(p(q));
            Xd_amp_post_m(i,j)=Xd_amp_post(p(q));
            Xd_pulse_prepulse_RMS_amp_m(i,j)=Xd_pulse_prepulse_RMS_amp(p(q));
            Xd_pulse_postpulse_RMS_amp_m(i,j)=Xd_pulse_postpulse_RMS_amp(p(q));
            post_rsquare_m(i,j)=post_rsquare(p(q));

        end
end


% Plot pulses / baselines / etc

u=size(Xd_pulse);
s=[logdata.data(1,iter_loc),a];

figure(1);
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(Xd_pulse(:,I(p),J(p)));
             hold on;             
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Xd Pulse'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(Xd_pulse(:,I(p),J(p)));
             hold on;
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Xd Pulse'])
       
      
    end
end
 
end

figure(2);
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(Xd_pre(:,I(p),J(p)));
             hold on;             
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Xd BASELINE'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(Xd_pre(:,I(p),J(p)));
             hold on;
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Xd BASELINE'])
       
      
    end
end
 
end

figure(3);
 
if N_F<15
   
for i = 1:N_F
    for j=1:N_k
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(N_F+1-i)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(j)));
       
       subplot(N_F,N_k,N_k*(i-1)+j)
       hold on;
       plot(Xd_pulse_postpulse_center(:,I(p),J(p)));
             hold on;             
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Xd Pulse Post-Pulse'])
       
      
    end
end
 
else
 
for i = 1:N_k
    for j=1:N_F
        
       [I,J]=ind2sub(s,find(F_rand==Ord_F(j)));
       
       p=ind2sub(s,find(k_rand(sub2ind(s,I,J))==Ord_k(i)));
       
       subplot(N_k,N_F,j+(i-1)*N_F)
       hold on;
       plot(Xd_pulse_postpulse_center(:,I(p),J(p)));
             hold on;
       title(['F=',num2str(F_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), ', k=',num2str(k_rand(sub2ind([u(2) u(3)],I(p),J(p))),3), 'Xd Pulse Post-Pulse'])
       
      
    end
end
 
end

choice2=menu('Force or Stiffness?','Force','Stiffness');
switch choice2
    case 1
        disp('Selected FORCE');
        figure(4);
        plot(Ord_F,Xd_freq_pre_m(:,1)); title('Frequency Pre-Pulse'); xlabel('Force');
        figure(5);
        plot(Ord_F,Xd_freq_post_m(:,1)); title('Frequency Post-Pulse'); xlabel('Force');
        figure(6);
        plot(Ord_F,Xd_amp_pre_m(:,1)); title('Amplitude Pre-Pulse'); xlabel('Force');        
        figure(7);
        plot(Ord_F,Xd_amp_post_m(:,1)); title('Amplitude Post-Pulse'); xlabel('Force');
        figure(8);
        plot(Ord_F,Xd_pulse_prepulse_RMS_amp_m(:,1)); title('RMS Amplitude Pre-Pulse'); xlabel('Force');     
        figure(9);
        plot(Ord_F,Xd_pulse_postpulse_RMS_amp_m(:,1)); title('RMS Amplitude Post-Pulse'); xlabel('Force'); 
    case 2
        disp('Selected STIFFNESS');
        figure(4);
        plot(Ord_k,Xd_freq_pre_m(1,:)); title('Frequency Pre-Pulse'); xlabel('Stiffness');
        figure(5);
        plot(Ord_k,Xd_freq_post_m(1,:)); title('Frequency Post-Pulse'); xlabel('Stiffness');
        figure(6);
        plot(Ord_k,Xd_amp_pre_m(1,:)); title('Amplitude Pre-Pulse'); xlabel('Stiffness');     
        figure(7);
        plot(Ord_k,Xd_amp_post_m(1,:)); title('Amplitude Post-Pulse'); xlabel('Stiffness');
        figure(8);
        plot(Ord_k,Xd_pulse_prepulse_RMS_amp_m(1,:)); title('RMS Amplitude Pre-Pulse'); xlabel('Stiffness');     
        figure(9);
        plot(Ord_k,Xd_pulse_postpulse_RMS_amp_m(1,:)); title('RMS Amplitude Post-Pulse'); xlabel('Stiffness'); 
end
    
else
    disp('Incorrect Choice');
end

%% Case 3 - Discrete Frequencies

if choice == 3
if choice2 ~= 4     
    
N=1;
 
min_freq = 5; % Fs/(length(Xd_pulse(:,1,1))/N);
max_freq = 100; %150;



  for j = 1:a
  for i = 1:(logdata.data(j,iter_loc))
 for p = 1:numfreq
stim_z(:,p,i,j)=smooth(stim(:,p,i,j),length(stim(:,p,i,j))/N);
stim_center(:,p,i,j) = stim(:,p,i,j) - stim_z(:,p,i,j);
Delta_in_z(:,p,i,j)=smooth(Delta_in(:,p,i,j),length(Delta_in(:,p,i,j))/N);
Delta_in_center(:,p,i,j) = Delta_in(:,p,i,j) - Delta_in_z(:,p,i,j);
%Xd_pre_z(:,i,j)=smooth(Xd_pre(:,i,j),length(Xd_pre(:,i,j))/N);
%Xd_post_z(:,i,j)=smooth(Xd_post(:,i,j),length(Xd_post(:,i,j))/N);
end
  end
  end
 

    clear q_stim psd_q1_stim psd_q2_stim

    
    T=1/Fs;
    NFFT = 2e4;
    
% STIMULUS, X
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
h_stim=psd(spectrum.periodogram('blackman'),stim_center(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',1e4);

L = length(stim_center(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j));
hfft  = fft(stim_center(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);

clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE


q = find((freq_stim(p)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));
q3 = findnearest(qfft_stim(:,1),freq_stim(p));

stim_amp(p,i,j)=sqrt(psd_q2_stim(q2,p,i,j)*psd_q1_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
stim_pow(p,i,j)=psd_q2_stim(q2,p,i,j);

stim_fft(p,i,j)=qfft_stim(q3,2);


end
end
end


    clear q_stim psd_q1Delta_stim psd_q2Delta_stim

% FIBER BASE, DELTA
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
h_stim=psd(spectrum.periodogram('blackman'),Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',2e4);

L = length(Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j));
hfft  = fft(Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_Delta(:,2)=hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);

clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1Delta_stim),:)=0;
end

psd_q1Delta_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2Delta_stim(:,p,i,j)=q_stim(:,2);    %psd PRE


q = find((freq_stim(p)-1)<psd_q1Delta_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1Delta_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2Delta_stim(q,p,i,j)==max(psd_q2Delta_stim(q,p,i,j))));
q3 = findnearest(qfft_stim(:,1),freq_stim(p));

Delta_in_amp(p,i,j)=sqrt(psd_q2Delta_stim(q2,p,i,j)*psd_q1Delta_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
Delta_in_pow(p,i,j)=psd_q2Delta_stim(q2,p,i,j);

Delta_in_fft(p,i,j)=qfft_Delta(q3,2);


end
end
end



% Calculate Sensitivity

% INPUT
ksf = 450e-6;       % Fiber stiffness (N/m)
gsf = 200e-9;       % Fiber damping (Ns/m)
amp = 5e-9;         % stimulus amplitude (m)

% Force = (ksf + i*w*gsf)*amp

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
sensitivity(p,i,j)=stim_amp(p,i,j)/(sqrt(ksf^2+(gsf*freq_stim(p)*2*pi)^2)*amp)*1e-9;
force_in(p,i,j) = (ksf + sqrt(-1)*freq_stim(p)*gsf)*Delta_in_fft(p,i,j);
res_func(p,i,j) = stim_fft(p,i,j)/force_in(p,i,j);
sensitivity_resfunc(p,i,j) = norm(res_func(p,i,j));
phase_resfunc(p,i,j) = -angle(res_func(p,i,j));
end
end
end

% Find means, STD, SEM --- Can repeat this for any other variable

for i = 1:numfreq
    for j = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
        mean_stim(i,j)=mean(stim_amp(i,j,:));
        std_stim(i,j)=std(stim_amp(i,j,:));
        sem_stim(i,j)=std_stim(i,j)/sqrt(a);
    end
end

% Find Q50
for j  = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for i= 2:numfreq-1

freq_stim_interp = min(freq_stim(1:numfreq)):0.01:max(freq_stim(1:numfreq));
mean_stim_interp = interp1(freq_stim(1:numfreq),mean_stim(:,j),freq_stim_interp);


amp_peak(i,j) = mean_stim(i,j);
amp_50(i,j) = amp_peak(i,j)/2;             % Q10dB... 10dB = 20*log10(A1/A2)
interp_left = mean_stim_interp(find(freq_stim_interp < freq_stim(i)));
interp_right = mean_stim_interp(find(freq_stim_interp > freq_stim(i)));

index1(i,j) = findnearest(interp_left,amp_50(i,j));
index2(i,j) = findnearest(interp_right,amp_50(i,j)) + length(interp_left);

Q50(i,j) = ((freq_stim_interp(index2(i,j))+freq_stim_interp(index1(i,j)))/2)/(freq_stim_interp(index2(i,j)) - freq_stim_interp(index1(i,j)));


end
end


figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),stim_amp(:,i,j),C{i});
title('Amplitude');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),sensitivity(:,i,j),C{i});
title('Sensitivity');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),sensitivity_resfunc(:,i,j),C{i});
title('Sensitivity from Response Function');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),phase_resfunc(:,i,j),C{i});
title('Phase from Response Function');
end
end

% Degree of Entrainment
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    stim_RMS(p,i,j)=std(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j));
    entrainment(p,i,j)=(stim_amp(p,i,j)*0.7)/stim_RMS(p,i,j);
end
end
end
%{
for j = 1:a
    if max(max(max(entrainment)))> 1
    entrainment(:,:,j)=entrainment(:,:,j)./max(max(max(entrainment(:,:,j))));
    end
end
%}
figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),entrainment(:,i,j),C{i});
title('Entrainment')
end
end


% Calculate Phase Lag (degrees)
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
[Cxy,f] = mscohere(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),[],[],[],Fs);    % Spectral Coherence.
Pxy     = cpsd(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),[],[],[],Fs);    % Fourier transform of cross-correlation.. Cross Power Spectral Density
phase   = -angle(Pxy)/pi*180;
q = findnearest(freq_stim(p),f);
%q=find(Cxy(f>freq_stim(p)/2&f<freq_stim(p)*2));
%q2=find(Cxy==max(Cxy(q)));
%q3=find(Pxy==max(Pxy(q)));
phase_lag(p,i,j)=phase(q);

coherence(p,i,j)=abs(Cxy(q));
cpsd2(p,i,j)=Pxy(q);

end
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),phase_lag(:,i,j),C{i});
title('Phase Lag (degrees), Cross Power Spectrum');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),coherence(:,i,j),C{i});
title('Coherence');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),-cpsd2(:,i,j),C{i});
title('Cross-Correlation Power');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),angle(cpsd2(:,i,j)),C{i});
title('Cross-Correlation Angle');
end
end


%{
% Calculate Time Lag
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
[c,lags]=xcorr(input(1:round(cycles/freq_stim(1)*1e4),p,i,j),stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),'coeff');
time_lag(p,i,j)=max(lags(c==max(c)))*1e-3;
phase_lag2(p,i,j)=time_lag(p,i,j)*freq_stim(p)*2*pi;
end
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim,phase_lag2(:,i,j),C{i});
title('Phase Lag (radians), Full Cross-Correlation');
end
end

% Find Transfer Function
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
[Txy(:,p,i,j),f2(:,p,i,j)]=tfestimate(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),[],[],[],Fs);
q = findnearest(f2,freq_stim(p));
phase_lag3(p,i,j)=angle(Txy(Txy==max(Txy(q))))*-1;
end
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim,phase_lag3(:,i,j),C{i});
title('Phase Lag (radians), Transfer Function');
end
end
%}
%{
% Reorder matrices

Ord_F=unique(F_rand);
Ord_k=unique(k_rand);

k=size(stim_amp);

stim_amp_sh=reshape(stim_amp,numfreq,k(2)*k(3));
stim_pow_sh=reshape(stim_pow,numfreq,k(2)*k(3));


for i=1:N_F
        %p=find(F_rand==Ord_F(N_F+1-i));
        p=find(F_rand==Ord_F(i));
        for j=1:N_k
            q=find(k_rand(p)==Ord_k(j));
 %           stim_amp_m(:,i,j)=stim_amp(:,p(q));
  %          stim_pow_m(:,i,j)=stim_pow(:,p(q));

        end
end
%}

else                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        
 
 disp('Amplitude Array');                                                       
                                                        
N=1;
 
min_freq = 5; % Fs/(length(Xd_pulse(:,1,1))/N);
max_freq = 100; %150;



  for j = 1:a
  for i = 1:(logdata.data(j,iter_loc))
 for p = 1:numfreq
stim_z(:,p,i,j)=smooth(stim(:,p,i,j),length(stim(:,p,i,j))/N);
stim_center(:,p,i,j) = stim(:,p,i,j) - stim_z(:,p,i,j);
Delta_in_z(:,p,i,j)=smooth(Delta_in(:,p,i,j),length(Delta_in(:,p,i,j))/N);
Delta_in_center(:,p,i,j) = Delta_in(:,p,i,j) - Delta_in_z(:,p,i,j);
%Xd_pre_z(:,i,j)=smooth(Xd_pre(:,i,j),length(Xd_pre(:,i,j))/N);
%Xd_post_z(:,i,j)=smooth(Xd_post(:,i,j),length(Xd_post(:,i,j))/N);
end
  end
  end
 

    clear q_stim psd_q1_stim psd_q2_stim

    
    T=1/Fs;
    NFFT = 2e4;
    
% STIMULUS, X
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
h_stim=psd(spectrum.periodogram('rectangular'),stim_center(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',1e4);

L = length(stim_center(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j));
hfft  = fft(stim_center(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);

clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE


q = find((freq_stim(p)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));
q3 = findnearest(qfft_stim(:,1),freq_stim(p));

stim_amp(p,i,j)=sqrt(psd_q2_stim(q2,p,i,j)*psd_q1_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
stim_pow(p,i,j)=psd_q2_stim(q2,p,i,j);

stim_fft(p,i,j)=qfft_stim(q3,2);


end
end
end


    clear q_stim psd_q1Delta_stim psd_q2Delta_stim

% FIBER BASE, DELTA
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
h_stim=psd(spectrum.periodogram('rectangular'),Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',2e4);

L = length(Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j));
hfft  = fft(Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_Delta(:,2)=hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);

clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1Delta_stim),:)=0;
end

psd_q1Delta_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2Delta_stim(:,p,i,j)=q_stim(:,2);    %psd PRE


q = find((freq_stim-1)<psd_q1Delta_stim(:,p,i,j)&1.05*(freq_stim+1)>psd_q1Delta_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2Delta_stim(q,p,i,j)==max(psd_q2Delta_stim(q,p,i,j))));
q3 = findnearest(qfft_stim(:,1),freq_stim);

Delta_in_amp(p,i,j)=sqrt(psd_q2Delta_stim(q2,p,i,j)*psd_q1Delta_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
Delta_in_pow(p,i,j)=psd_q2Delta_stim(q2,p,i,j);

Delta_in_fft(p,i,j)=qfft_Delta(q3,2);


end
end
end



% Calculate Sensitivity

% INPUT
ksf = 450e-6;       % Fiber stiffness (N/m)
gsf = 200e-9;       % Fiber damping (Ns/m)
amp = 5e-9;         % stimulus amplitude (m)

% Force = (ksf + i*w*gsf)*amp

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
sensitivity(p,i,j)=stim_amp(p,i,j)/(sqrt(ksf^2+(gsf*freq_stim*2*pi)^2)*amp)*1e-9;
force_in(p,i,j) = (ksf + sqrt(-1)*freq_stim*gsf)*Delta_in_fft(p,i,j);
res_func(p,i,j) = stim_fft(p,i,j)/force_in(p,i,j);
sensitivity_resfunc(p,i,j) = norm(res_func(p,i,j));
phase_resfunc(p,i,j) = -angle(res_func(p,i,j));
end
end
end

% Find means, STD, SEM --- Can repeat this for any other variable

for i = 1:numfreq
    for j = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
        mean_stim(i,j)=mean(stim_amp(i,j,:));
        std_stim(i,j)=std(stim_amp(i,j,:));
        sem_stim(i,j)=std_stim(i,j)/sqrt(a);
    end
end

%{
% Find Q50
for j  = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for i= 1:numfreq

amp_stim_interp = min(amp_stim(1:numfreq)):0.01:max(freq_stim);
mean_stim_interp = interp1(amp_stim(1:numfreq),mean_stim(:,j),amp_stim_interp);


amp_peak(i,j) = mean_stim(i,j);
amp_50(i,j) = amp_peak(i,j)/2;             % Q10dB... 10dB = 20*log10(A1/A2)
interp_left = mean_stim_interp(find(amp_stim_interp < amp_stim(i)));
interp_right = mean_stim_interp(find(amp_stim_interp > amp_stim(i)));

index1(i,j) = findnearest(interp_left,amp_50(i,j));
index2(i,j) = findnearest(interp_right,amp_50(i,j)) + length(interp_left);

Q50(i,j) = ((amp_stim_interp(index2(i,j))+freq_stim_interp(index1(i,j)))/2)/(amp_stim_interp(index2(i,j)) - amp_stim_interp(index1(i,j)));


end
end

%}
figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; loglog(amp_stim(1:numfreq),stim_amp(:,i,j),C{i});
title('Amplitude');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; loglog(amp_stim(1:numfreq),sensitivity(:,i,j),C{i});
title('Sensitivity');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; loglog(amp_stim(1:numfreq),sensitivity_resfunc(:,i,j),C{i});
title('Sensitivity from Response Function');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; loglog(amp_stim(1:numfreq),phase_resfunc(:,i,j),C{i});
title('Phase from Response Function');
end
end

% Degree of Entrainment
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    stim_RMS(p,i,j)=std(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j));
    entrainment(p,i,j)=(stim_amp(p,i,j)*0.7)/stim_RMS(p,i,j);
end
end
end
%{
for j = 1:a
    if max(max(max(entrainment)))> 1
    entrainment(:,:,j)=entrainment(:,:,j)./max(max(max(entrainment(:,:,j))));
    end
end
%}
figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(amp_stim,entrainment(:,i,j),C{i});
title('Entrainment')
end
end

% Vector Strength
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    chilb = hilbert(input(1:stim_time(p),p,i,j));
    cdelt = angle(chilb);
    xhilb = hilbert(stim_center(1:stim_time(p),p,i,j));
    xdelt = angle(xhilb);
    
    sphase = cdelt - xdelt;
    
    vec_length(p,i,j) = circ_r(sphase);
    
    figure(2*a+j)
    subaxis(logdata.data(j,iter_loc),numfreq,p+(i-1)*numfreq,'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
    set(gca, 'LooseInset', get(gca,'TightInset'))
    polar(0, 0.1); hold on;
    circ_plot(sphase,'hist',[],N,true,true,'linewidth',2,'color','r');
end
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(amp_stim,vec_length(:,i,j),C{i});
title('Vector Strength')
end
end
%{
% Calculate Phase Lag (degrees)
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
[Cxy,f] = mscohere(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),[],[],[],Fs);    % Spectral Coherence.
Pxy     = cpsd(stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),[],[],[],Fs);    % Fourier transform of cross-correlation.. Cross Power Spectral Density
phase   = -angle(Pxy)/pi*180;
q = findnearest(freq_stim,f);
%q=find(Cxy(f>freq_stim/2&f<freq_stim*2));
%q2=find(Cxy==max(Cxy(q)));
%q3=find(Pxy==max(Pxy(q)));
phase_lag(p,i,j)=phase(q);

coherence(p,i,j)=abs(Cxy(q));
cpsd2(p,i,j)=Pxy(q);

end
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),phase_lag(:,i,j),C{i});
title('Phase Lag (degrees), Cross Power Spectrum');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),coherence(:,i,j),C{i});
title('Coherence');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),-cpsd2(:,i,j),C{i});
title('Cross-Correlation Power');
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim(1:numfreq),angle(cpsd2(:,i,j)),C{i});
title('Cross-Correlation Angle');
end
end

%}
%{
% Calculate Time Lag
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
[c,lags]=xcorr(input(1:round(cycles/freq_stim(1)*1e4),p,i,j),stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),'coeff');
time_lag(p,i,j)=max(lags(c==max(c)))*1e-3;
phase_lag2(p,i,j)=time_lag(p,i,j)*freq_stim*2*pi;
end
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim,phase_lag2(:,i,j),C{i});
title('Phase Lag (radians), Full Cross-Correlation');
end
end

% Find Transfer Function
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
[Txy(:,p,i,j),f2(:,p,i,j)]=tfestimate(Delta_in_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),stim_center(1:round(cycles/freq_stim(1)*1e4),p,i,j),[],[],[],Fs);
q = findnearest(f2,freq_stim);
phase_lag3(p,i,j)=angle(Txy(Txy==max(Txy(q))))*-1;
end
end
end

figure();
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
C = {'k','b','r','g','y','k+-','b+-'};
subplot(a,1,j);
hold on; plot(freq_stim,phase_lag3(:,i,j),C{i});
title('Phase Lag (radians), Transfer Function');
end
end
%}
%{
% Reorder matrices

Ord_F=unique(F_rand);
Ord_k=unique(k_rand);

k=size(stim_amp);

stim_amp_sh=reshape(stim_amp,numfreq,k(2)*k(3));
stim_pow_sh=reshape(stim_pow,numfreq,k(2)*k(3));


for i=1:N_F
        %p=find(F_rand==Ord_F(N_F+1-i));
        p=find(F_rand==Ord_F(i));
        for j=1:N_k
            q=find(k_rand(p)==Ord_k(j));
 %           stim_amp_m(:,i,j)=stim_amp(:,p(q));
  %          stim_pow_m(:,i,j)=stim_pow(:,p(q));

        end
end
%}   
  

end
else
    disp('Incorrect Choice');
end

%% Case 4 - Frequency Sweep

if choice == 4
    
    
    
    
    
    
else
    disp('Incorrect Choice');
end

%% Case 5 - Ramp

if choice == 5
    
    
    
    
    
    
else
    disp('Incorrect Choice');
end
