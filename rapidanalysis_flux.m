%%
figure;
%plot(XdInput1Displacementnm11+50);
hold on;
plot(XdInput1Displacementnm12+50);
plot(XdInput1Displacementnm13+0);
%plot(XdInput1Displacementnm14-100);
plot(XdInput1Displacementnm15-50);
plot(XdInput1Displacementnm16-100);
plot(XdInput1Displacementnm17-150);
%plot(XdInput1Displacementnm18-300);
plot(XdInput1Displacementnm19-200);

%%
%close all;
N2 = 1000;
N= 5;
Fs=5e3;
start=2e4;
endl=4.5e4;

[pxx1 f] = pmtm(XdInput1Displacementnm11(start:endl)-smooth(XdInput1Displacementnm11(start:endl),N2),N,[],Fs);
[pxx2 f] = pmtm(XdInput1Displacementnm12(start:endl)-smooth(XdInput1Displacementnm12(start:endl),N2),N,[],Fs);
[pxx3 f] = pmtm(XdInput1Displacementnm13(start:endl)-smooth(XdInput1Displacementnm13(start:endl),N2),N,[],Fs);
[pxx4 f] = pmtm(XdInput1Displacementnm14(start:endl)-smooth(XdInput1Displacementnm14(start:endl),N2),N,[],Fs);
[pxx5 f] = pmtm(XdInput1Displacementnm15(start:endl)-smooth(XdInput1Displacementnm15(start:endl),N2),N,[],Fs);
[pxx6 f] = pmtm(XdInput1Displacementnm16(start:endl)-smooth(XdInput1Displacementnm16(start:endl),N2),N,[],Fs);
[pxx7 f] = pmtm(XdInput1Displacementnm17(start:endl)-smooth(XdInput1Displacementnm17(start:endl),N2),N,[],Fs);
[pxx8 f] = pmtm(XdInput1Displacementnm18(start:endl)-smooth(XdInput1Displacementnm18(start:endl),N2),N,[],Fs);
[pxx9 f] = pmtm(XdInput1Displacementnm19(start:endl)-smooth(XdInput1Displacementnm19(start:endl),N2),N,[],Fs);

figure;
C = {'k','b','g','r','c','m','y','b--','g--'};
%loglog(f,pxx1,C{1});hold on;
loglog(f,pxx2,C{2});hold on;
loglog(f,pxx3,C{3});hold on;
%loglog(f,pxx4,C{4});hold on;
loglog(f,pxx5,C{5});hold on;
loglog(f,pxx6,C{6});hold on;
loglog(f,pxx7,C{7});hold on;
%loglog(f,pxx8,C{8});hold on;
loglog(f,pxx9,C{9});hold on;

%% RMS Amplitude

N2 = 1000;
N= 5;
Fs=5e3;
start=0.5e4;
endl=2.5e4;

RMSamp(1) = std(XdInput1Displacementnm11(start:endl)-smooth(XdInput1Displacementnm11(start:endl),N2));
RMSamp(2) = std(XdInput1Displacementnm12(start:endl)-smooth(XdInput1Displacementnm12(start:endl),N2));
RMSamp(3) = std(XdInput1Displacementnm13(start:endl)-smooth(XdInput1Displacementnm13(start:endl),N2));
RMSamp(4) = std(XdInput1Displacementnm14(start:endl)-smooth(XdInput1Displacementnm14(start:endl),N2));
RMSamp(5) = std(XdInput1Displacementnm15(start:endl)-smooth(XdInput1Displacementnm15(start:endl),N2));
RMSamp(6) = std(XdInput1Displacementnm16(start:endl)-smooth(XdInput1Displacementnm16(start:endl),N2));
RMSamp(7) = std(XdInput1Displacementnm17(start:endl)-smooth(XdInput1Displacementnm17(start:endl),N2));
RMSamp(8) = std(XdInput1Displacementnm18(start:endl)-smooth(XdInput1Displacementnm18(start:endl),N2));
RMSamp(9) = std(XdInput1Displacementnm19(start:endl)-smooth(XdInput1Displacementnm19(start:endl),N2));
