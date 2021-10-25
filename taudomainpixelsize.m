function [tau] = taudomainpixelsize(numpixels,sigma)
% González-Siu L.O., 13/03/2020
% Stokes polarimetry by channeled polarimetry - Simulation

% This function returns a (tau) vector containing the corresponding
% positions in the tau domain (Inverse Fourier Space).
% Given the number of pixels and the range [sigma(1),sigma(end)] in [m]

tau = (-numpixels/2:numpixels/2-1).*(1/(sigma(end)-sigma(1)))*1e6;        % [um]

end