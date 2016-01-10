%% simulation parameters
nTimepoints=50000; % total number of time points
nOccSig1=100; % number of occurances of signal 1
nOccSig2=100; % number of occurances of signal 2
sigLength=1000; % signals length in timepoints
noiseAmplitude=10^-3; % random walk noise magnitude (try 10^-3)
baseline=7;

%% generate 'true signals'
trueSig1=cos(linspace(-pi,pi,sigLength))'*.5+.5;
trueSig2=[zeros(sigLength*.3,1);ones(sigLength*.4,1);zeros(sigLength*.3,1)];

% choose random event timings
occSig1=sort(ceil(rand(nOccSig1,1)*(nTimepoints-sigLength)));
occSig2=sort(ceil(rand(nOccSig2,1)*(nTimepoints-sigLength)));

%% generate observed timecourse
% random walk noise
noise=cumsum(randn(nTimepoints,1)).*noiseAmplitude+baseline; 
observedTimecourse=noise;

% add signal 1
for iOcc=1:nOccSig1
    curOcc=occSig1(iOcc); % current event occurance
    observedTimecourse(curOcc:(curOcc+sigLength-1))=observedTimecourse(curOcc:(curOcc+sigLength-1))+trueSig1;
end

% add signal 2
for iOcc=1:nOccSig2
    curOcc=occSig2(iOcc); % current event occurance
    observedTimecourse(curOcc:(curOcc+sigLength-1))=observedTimecourse(curOcc:(curOcc+sigLength-1))+trueSig2;
end

%% standard erp analysis
% estimate signal 1
M=nan(sigLength,nOccSig1);
for iOcc=1:nOccSig1
    curSeg=occSig1(iOcc):(occSig1(iOcc)+sigLength-1);
    M(:,iOcc)=observedTimecourse(curSeg);
end
erpEstimatedSig1=mean(M,2)-mean(observedTimecourse); clear M

% estimate signal 2
M=nan(sigLength,nOccSig2);
for iOcc=1:nOccSig2
    curSeg=occSig2(iOcc):(occSig2(iOcc)+sigLength-1);
    M(:,iOcc)=observedTimecourse(curSeg);
end
erpEstimatedSig2=mean(M,2)-mean(observedTimecourse); clear M

%% overlap correction 

latencies=[occSig1;occSig2];
types=[1*ones(numel(occSig1),1);2*ones(numel(occSig1),1)];
glm_estimates = rerp( observedTimecourse, latencies,types,sigLength);
glmEstimatedSig1=glm_estimates(:,1);
glmEstimatedSig2=glm_estimates(:,2);


%% plot results
figure
subplot(2,2,1);
plot(erpEstimatedSig1,'g','linewidth',2); ylim([-1 2]); hold all;
plot(trueSig1,'--k'); hold off;
title ('signal1 - erp');

subplot(2,2,3);
plot(erpEstimatedSig2,'g','linewidth',2); ylim([-1 2]); hold all;
plot(trueSig2,'--k'); hold off;
title ('signal2 - erp');


subplot(2,2,2);
plot(glmEstimatedSig1,'r','linewidth',2); ylim([-1 2]); hold all;
plot(trueSig1,'--k'); hold off;
title ('signal1 - glm');

subplot(2,2,4);
plot(glmEstimatedSig2,'r','linewidth',2); ylim([-1 2]); hold all;
plot(trueSig2,'--k'); hold off;
title ('signal2 - glm');

