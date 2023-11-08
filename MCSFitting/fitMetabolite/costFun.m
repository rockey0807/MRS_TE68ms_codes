function F = costFun(coeff,xdata,ydata,basisFID,nBasis,mrsPars)

fidFit = zeros(1,mrsPars.np);

for i = 1 : nBasis
    
    tmp = coeff(i) * basisFID(i,:);   
    fidFit = fidFit + tmp;
    
    clear tmp
    
end

%***************************************************
% Fourier transformation (optional zero filling)
%***************************************************
specFit = fftshift(fft(fidFit,mrsPars.zf));

specFit = reshape(specFit,length(specFit),1);
specFitROI = real(specFit(xdata)); 

F = specFitROI;


