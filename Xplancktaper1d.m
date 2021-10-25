function [fPT] = Xplancktaper1d(W,eps)
% Function to define a 1-dimensional Planck-Taper filter, or bump function,
% As defined in
% https://en.wikipedia.org/wiki/Window_function#Planck-taper_window
% W has to be an even number

if mod(W,2) == 1
  W = W + 1;
end

xarray = 1:W;

fPT = zeros(1,length(xarray));
n = 1:floor(eps*W);
fPT(n) = 1./(1 + exp(eps*W./n - eps*W./(eps*W-n)));
fPT(ceil(eps*W):ceil(W/2)) = 1;
fPT = fPT + flip(fPT);

end