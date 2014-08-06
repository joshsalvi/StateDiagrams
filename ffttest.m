tic

clear H_fft P_fft CI_fft stats_fft fft_sd p freq_comp amp_comp freq_test

window_size = 17;
overlap = 0.5;

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
    for k = 1:floor(length(psd_q1fft)/(window_size))-1+floor(length(psd_q1fft)/(window_size))*(1-overlap)
%        smoothed_psd(:,i,j) = smooth(abs(psd_q2fft(:,i,j)),window_size);
if k > 1
        [H_fft(i,j,k),P_fft(i,j,k),CI_fft(:,i,j,k),stats_fft(:,i,j,k)]=ttest(abs(psd_q2fft((k-1)*round(window_size*overlap)-round(window_size*overlap)+1:(k)*round(window_size*overlap)-round(window_size*overlap),i,j)),abs(psd_q2fft((k)*round(window_size*overlap)-round(window_size*overlap)+1:(k+1)*round(window_size*overlap)-round(window_size*overlap),i,j)));
else
        [H_fft(i,j,k),P_fft(i,j,k),CI_fft(:,i,j,k),stats_fft(:,i,j,k)]=ttest(abs(psd_q2fft((k-1)*round(window_size*overlap)+1:(k)*round(window_size*overlap),i,j)),abs(psd_q2fft((k)*round(window_size*overlap)-round(window_size*overlap)+1:(k+1)*round(window_size*overlap)-round(window_size*overlap),i,j)));
end    
        fft_sd(i,j,k)=stats_fft(:,i,j,k).sd;

    end
end
end

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
    p = find(H_fft(i,j,:)==1);
    for k = 1:length(p)
    freq_comp(:,i,j,k)=psd_q1fft(p(k)*round(window_size*overlap)-round(window_size*overlap)+1:(p(k)+1)*round(window_size*overlap)-round(window_size*overlap),i,j);
    amp_comp(:,i,j,k)=psd_q2fft(p(k)*round(window_size*overlap)-round(window_size*overlap)+1:(p(k)+1)*round(window_size*overlap)-round(window_size*overlap),i,j);
    end    
end
end

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
    freq_test(i,j)=0;
    for k = 1:size(freq_comp,4)
        if isempty(find(round(freq_comp(:,i,j,k))==round(Xd_freq(i,j))))==0
            freq_test(i,j)=1;
        end
    end
end
end


toc