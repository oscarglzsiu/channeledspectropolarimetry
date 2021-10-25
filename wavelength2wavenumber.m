function [sigma,I_sigma,deltasigma] = wavelength2wavenumber(lambda,I_lambda,N_sigma)
% by Gonzalez-Siu, Luis Oscar (12/08/2021)
% This function converts a spectrum recorded in wavelength (lambda in [nm])
% to a uniformly spaced spectrum in wavenumber (sigma in [m-1]).
sigma_nonuniform = 1e9./lambda;
sigmamin = sigma_nonuniform(end);
sigmamax = sigma_nonuniform(1);
deltasigma = (sigmamax - sigmamin) / N_sigma;
% Uniformly spaced wavenumber vector
sigma = sigmamin:deltasigma:sigmamax-deltasigma;
% Interpolation of a uniformly spaced spectrum in wavenumber
I_sigma = interp1(flip(sigma_nonuniform),flip(I_lambda),sigma);
end