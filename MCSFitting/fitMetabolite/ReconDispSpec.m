function ReconDispSpec(lcmFit,spec,basisFID,mrsPars,lcmPars,nBasis)

%
fidFit = zeros(1,mrsPars.np);
fidFitInd = zeros(nBasis,mrsPars.np);

for i = 1 : nBasis
    
    fidFitInd(i,:) = lcmFit(i) * basisFID(i,:);   
    fidFit = fidFit + fidFitInd(i,:);
    
end

%***************************************************
% Fourier transformation (optional zero filling)
%***************************************************
specFit = fftshift(fft(fidFit,mrsPars.zf));

specFitInd = zeros(nBasis,mrsPars.zf);
for i = 1 : nBasis
    
    specFitInd(i,:) = fftshift(fft(fidFitInd(i,:),mrsPars.zf));
    
end

[~,lowIdx] = min(abs(mrsPars.ppm - lcmPars.fitRange(1)));
[~,highIdx] = min(abs(mrsPars.ppm - lcmPars.fitRange(2)));

xdata = lowIdx : highIdx;

spec = reshape(spec,1,length(spec));

rawSpec = real(spec(xdata)); 
fitSpec = real(specFit(xdata)); 

residual = rawSpec - fitSpec;

ppm = mrsPars.ppm(xdata);

% Display and save fitting results
fontsize = 11;
LW = 1.25;

minY = min(rawSpec);
minY = minY - 50;
maxY = max(residual) + max(rawSpec) + 70;

figure('units','inch','Position',[2 2 7 2.5]);

subplot(1,2,1); plot(ppm,real(rawSpec),'LineWidth',LW,'Color',[0 0 0])
hold on; plot(ppm,real(fitSpec),'LineWidth',LW,'LineStyle','--','Color',[1 0 0])
plot(ppm,residual+max(rawSpec)+50,'LineWidth',LW,'Color',[0 0 0])

xlim([lcmPars.fitRange(1) lcmPars.fitRange(2)])
xticks([2.0 2.5 3.0 3.5 4.0]); xticklabels({'2.0','2.5','3.0','3.5','4.0'})
set(gca,'XDir','reverse','FontSize',fontsize,'LineWidth',LW,'box','off','ycolor','none')
ylim([minY maxY])

%
subplot(1,2,2); plot(ppm,real(rawSpec),'LineWidth',LW,'Color',[0 0 0])

for i = 1 : nBasis
    tmp = real(specFitInd(i,xdata));
    
    hold on; plot(ppm,tmp,'LineWidth',LW)
    clear tmp
end

xlim([lcmPars.fitRange(1) lcmPars.fitRange(2)])
xticks([2.0 2.5 3.0 3.5 4.0]); xticklabels({'2.0','2.5','3.0','3.5','4.0'})
set(gca,'XDir','reverse','FontSize',fontsize,'LineWidth',LW,'box','off','ycolor','none')
ylim([minY maxY])
