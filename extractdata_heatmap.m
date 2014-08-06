% Extract data from FIG file
clear all; close all;

open('/Users/joshsalvi/Documents/Lab/Lab/Original/Paper/Raw Data/State Space Analysis/SmallCellResults/FreqFishmedian0.67sDetrend3smin1Hzmin.fig');
  
D=get(gca,'Children'); %get the handle of the line object
XData=get(D,'XData'); %get the x data
YData=get(D,'YData'); %get the y data
%ZData=get(D,'ZData');  %get z data
CData=get(D,'CData');   %get colormap data


%% FOR HEATMAPS

% Sum data along axes of heatmap
for i = 1:7
for j = 1:7
CDatasum1(i) = mean(CData(i,isnan(CData(i,:))==0));
CDatasum2(j) = mean(CData(isnan(CData(:,j))==0,j));
end
end
CDatasum1(isnan(CDatasum1)==1)=0;
CDatasum2(isnan(CDatasum2)==1)=0;

ampl=0;     % 0 for freq, 1 for ampl
freq=1;
if ampl ==1
figure(2);subplot(2,1,1);area(XData(1:7),CDatasum2,'FaceColor',[0.6 0 0]); xlabel('Stiffness (µN/m)'); ylabel('Mean Amplitude')
subplot(2,1,2);area(YData(1:7),CDatasum1,'FaceColor',[0.6 0 0]); xlabel('Force (pN)'); ylabel('Mean Amplitude');
else if freq ==1
figure(2);subplot(2,1,1);area(XData(1:7),CDatasum2,'FaceColor',[0 0 0.6]); xlabel('Stiffness (µN/m)'); ylabel('Mean Frequency')
subplot(2,1,2);area(YData(1:7),CDatasum1,'FaceColor',[0 0 0.6]); xlabel('Force (pN)'); ylabel('Mean Frequency'); 
    end
end

