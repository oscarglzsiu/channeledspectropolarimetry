function B = birefringenceSW(crystalid,lambda)
% by Gonzalez-Siu, Luis Oscar (12/08/2021)
% This function returns the birefringence profile B(lambda in nm) given the
% crystalid and the wavelenght array lambda in nm
switch crystalid
%  Case 1 and 2 use parameters and model taken from Ghosh (1999).
%  Ghosh's model considers wavelength in micrometers (um)
case 1
  % Quartz
  H = 0.78890253e-3;
  I = 8.04095323e-3;
  G = 1.37254429e-2;
  J = 10.1933186e-3;
  L = 64;
  lambdaum = lambda * 1e-3;
  B = H + I.*lambdaum.^2./(lambdaum.^2-G)+J.*lambdaum.^2./(lambdaum.^2-L);
case 2
  % Calcite
  H = -29.435688e-3;
  I = -134.804456e-3;
  G = 2.17641576e-2;
  J = -294.96110e-3;
  L = 80;
  lambdaum = lambda * 1e-3;
  B = H + I.*lambdaum.^2./(lambdaum.^2-G)+J.*lambdaum.^2./(lambdaum.^2-L);
end
end
