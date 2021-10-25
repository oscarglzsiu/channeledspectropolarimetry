function phi = retardance(birefringence,thickness,lambda)
% by Gonzalez-Siu, Luis Oscar (12/08/2021)
% This function returns the retardance phi(lambda) given the retarder's
% birefringence, thickness in m, and wavelength band lambda in nm
phi = 2*pi .* birefringence * thickness * 1e9 ./ lambda;
end