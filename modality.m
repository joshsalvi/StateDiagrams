bimod_gauss = fittype('p*1/(sigma1*sqrt(2*pi))*exp(-0.5*((x-mu1)/sigma1))^2+(p-1)*1/(sigma2*sqrt(2*pi))*exp(-0.5*((x-mu2)/sigma2))^2');

clear a b AIC obj k options model numComponents


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