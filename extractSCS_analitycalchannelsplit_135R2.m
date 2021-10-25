function [extractedS,extractedSn] = extractSCS_analitycalchannelsplit_135R2(...
  AChannels,A0index,Qmatrix,channelposition,...
  phi1,phi2,thicknessratio,omega)
% by Gonzalez-Siu, Luis Oscar (12/08/2021)
% This function returns the extracted Stokes parameters using the
% Analytical Channel Splitting method described by Alenin and Tyo (2012).
% Pseudo-invert the Q-matrix
Qinv = pinv(Qmatrix);
% Fourier Transform of the C-'vector' (C channels matrix) along dim = 2
% (i.e., along the rows) and reverse apodization.
fftCChanels = fftshift( fft( fftshift(AChannels,2) ,[],2),2) ./ omega;
% Phase cancelation ------------------------------------------------------
% Channels located at +/-tau_2
fftCChanels(A0index + channelposition(1),:) = ...
    fftCChanels(A0index + channelposition(1),:) .* exp(1i.*phi2);
fftCChanels(A0index - channelposition(1),:) = ...
    fftCChanels(A0index - channelposition(1),:) .* exp(-1i.*phi2);
% Channels located at +/-(tau_1+tau_2)
fftCChanels(A0index + channelposition(2),:) = ...
    fftCChanels(A0index + channelposition(2),:) .* exp(1i.*(phi1+phi2));
fftCChanels(A0index - channelposition(2),:) = ...
    fftCChanels(A0index - channelposition(2),:) .* exp(-1i.*(phi1+phi2));
% Channels located at +/-(tau_1-tau_2), only when configuration ~= (2,1)
if (thicknessratio(1) / thicknessratio(2) ~= 2)
fftCChanels(A0index + channelposition(3),:) = ...
    fftCChanels(A0index + channelposition(3),:) .* exp(1i.*(phi1-phi2));
fftCChanels(A0index - channelposition(3),:) = ...
    fftCChanels(A0index - channelposition(3),:) .* exp(-1i.*(phi1-phi2));
end
% Extract the Stokes parameters using the Q-matrix
extractedS = real(Qinv * fftCChanels);
% Normalization of Stokes parameters S1/S0, S2/S0, and S3/S0
extractedSn = extractedS(2:end,:)./extractedS(1,:);
end