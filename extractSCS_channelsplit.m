function [extractedS,extractedSn] = extractSCS_channelsplit...
  (AChannels,A0index,channelposition,...
  numpixels,phi1,phi2,thicknessratio,omega)
% by Gonzalez-Siu, Luis Oscar (12/08/2021)
% This function returns the extracted Stokes parameters using the Channel
% Splitting method described by Oka and Kato (1999).
extractedS = zeros(4,numpixels);
% Apply the Fourier Transform to the channel located at (tau_1+tau_2)
Fsum = fftshift( fft( fftshift(...
    AChannels(A0index + channelposition(2),:) )));
supportS23 = -8 * (Fsum .* exp(1i.*(phi1+phi2)));
extractedS(3,:) = real(supportS23);
extractedS(4,:) = imag(supportS23);
% Apply the Fourier Transform to the channel located at tau_2
F2 = fftshift( fft( fftshift(...
    AChannels(A0index + channelposition(1),:) )));
extractedS(2,:) = 4 * F2 .* exp(1i.*phi2);
if (thicknessratio(1) / thicknessratio(2) == 2)
    extractedS(2,:) = extractedS(2,:) - 0.5 * supportS23;
end
% Apply the Fourier Transform to the channel located at tau = 0
F0 = fftshift( fft( fftshift( AChannels(A0index,:) )));
extractedS(1,:) = 2*F0;
% Reverse apodization
extractedS = real(extractedS ./ omega);
% Normalization of Stokes parameters S1/S0, S2/S0, and S3/S0
extractedSn = extractedS(2:end,:) ./ extractedS(1,:);
end