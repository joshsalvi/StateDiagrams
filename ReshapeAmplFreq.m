% This creates the 1D vectors necessary for correlations between state
% diagrams across different cells.

%% RMS Magnitude

amplgrid(:,8)=[];
amplgrid(8,:)=[];

RMSampl = reshape(amplgrid,1,49);

save('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/Gentamicin/2014-08-05.01/Ear 1/Cell 7/RMSampl4sec2Hzmin.mat','RMSampl');


%% Amplitude and Frequency

amplgrid(:,8)=[];
amplgrid(8,:)=[];

freqgrid(:,8)=[];
freqgrid(8,:)=[];

ampl = reshape(amplgrid,1,49);
freq = reshape(freqgrid,1,49);

save('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/Controls/Gentamicin/2014-08-05.01/Ear 1/Cell 7/amplfreq3sec2Hzmin.mat','ampl','freq');