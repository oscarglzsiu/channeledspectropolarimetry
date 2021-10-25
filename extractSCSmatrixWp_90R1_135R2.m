function [extractedS,extractedSn] = extractSCSmatrixWp_90R1_135R2...
  (lambda,intensity,numpixels,numframes,thicknessvector,B1,B2)
% by Gonzalez-Siu, Luis Oscar (12/08/2021)
% This function extracts the Stokes parameters as a function of wavelength
% (lambda) using characteristic sub-matrices Wp representative of frames
% taken from the irradiance spectrum. With this method, we avoid the
% convertion from a conventional spectrum in wavelength to another
% non-uniform spectrum in wavenumber.
% B1 = birefringenceswitch(crystalid(1),lambda); % Birefringence
% B2 = birefringenceswitch(crystalid(2),lambda);
phi1 = retardance(B1,thicknessvector(1),lambda); % Retardance
phi2 = retardance(B2,thicknessvector(2),lambda);
W = mueller_SCSW_90R1_135R2(phi1,phi2); % Characteristic matrix of the SCS
ppframe = floor(numpixels/numframes); % Pixels per frame
extractedS = zeros(4,numframes); %Initialize the extracted Stokes vector
for i = 1:numframes
%     Take sections from the W-matrix for each frame
%     Wp = W((1:ppframe)+(i-1)*ppframe,:);
%     Pseudo-inverse of matrix Wp using Singular Value Decomposition (SVD) 
%     Wpinv = pinv( W((1:ppframe)+(i-1)*ppframe,:) );
%     Extract Stokes vector for i-th frame
    extractedS(:,i) = pinv( W((1:ppframe)+(i-1)*ppframe,:) ) *...
      intensity((1:ppframe)+(i-1)*ppframe).';
end
% Normalization of Stokes parameters S1/S0, S2/S0, and S3/S0
extractedSn = extractedS(2:end,:) ./ extractedS(1,:);
end