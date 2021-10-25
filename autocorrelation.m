function [functionC,magnitudeC] = autocorrelation(intensity)
% by Gonzalez-Siu, Luis Oscar (25/10/2019)
% This function returns the autocorrelation function C(tau) and its
% magnitude |C(tau)| of an intensity spectrum given by the inverse Fourier
% Transform of the intensity.
functionC = fftshift( ifft( fftshift( intensity )));
magnitudeC = abs(functionC);
end