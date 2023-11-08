clear; close all; clc

% MRS parameters
mrsPars.freq = 127.73;
mrsPars.sw = 5;
mrsPars.np = 4096;
mrsPars.rfoffset = 4.67;
mrsPars.basisRFoffset = 4.67;
mrsPars.zf = 2 * mrsPars.np;
mrsPars.ff = -0.5*mrsPars.sw : mrsPars.sw/(mrsPars.zf-1) : 0.5*mrsPars.sw;
mrsPars.ppm = mrsPars.rfoffset + 1000*mrsPars.ff/mrsPars.freq;

% LCM fitting parameters
lcmPars.fitRange = [1.8 4.2];
lcmPars.pointBsl = 10;
lcmPars.maxiter = 30;

% MCS parameters
numData = 20;

% Read data
dataLoad = '\Share\Data\MCS\';

% Load basis data
name = {'NAA','Glu','NAAG','tCr','tCho','mI','Asp','GABA',...
    'Glc','Gln','Gly','GSH','PE','sI','Tau'};

basisDir = '\Share\MCSshare\BasisSet\';
nBasis = length(name);

for i = 1 : nBasis
    
    fileName = [basisDir,name{i}];
    basisFID(i,:) = readData(fileName);
    
end

%***************************************
% set initial, low and upper boundaries
%***************************************
LB = zeros(1,nBasis);
UB = Inf * ones(1,nBasis);
init = ones(1,nBasis);

%********************************************
% Determine data points used for LCM fitting
%********************************************
[~,lowIdx] = min(abs(mrsPars.ppm - lcmPars.fitRange(1)));
[~,highIdx] = min(abs(mrsPars.ppm - lcmPars.fitRange(2)));

fitRange = lowIdx : highIdx;
xdata = reshape(fitRange,length(fitRange),1);

%*********************************************
% Set optimization parameters for LCM fitting
%*********************************************
tolfun = 1e-10;
maxfuneval = lcmPars.maxiter * (length(init)+1);

opt = optimset('Tolfun',tolfun,'MaxFunEval',maxfuneval,'LargeScale','on','DerivativeCheck','off','FinDiffType','central', ...
    'MaxIter',lcmPars.maxiter,'Display','iter','DiffMinChange',1e-8,'UseParallel','always');

%******************
% Start fitting
%******************
mcsAmp = zeros(numData,nBasis);

for i = 1 : numData
    
    FID = readData([dataLoad 'FID' num2str(i) '.dat']);
    spec = fftshift(fft(FID,mrsPars.zf));
    spec = reshape(spec,length(spec),1);
    
    ydata = real(spec(xdata));
    
    [lcmFit,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit('costFun',init,xdata,ydata,LB,UB,opt,ydata,basisFID,nBasis,mrsPars);
    
    if i < 10
        ReconDispSpec(lcmFit,spec,basisFID,mrsPars,lcmPars,nBasis);
        pause(1);
    end
    
    mcsAmp(i,:) = lcmFit;
    
    clear noiseComp noiseSpec ydata lcmFit
    
end

% save
save([dataLoad 'mcsAmp.mat'],'mcsAmp');
