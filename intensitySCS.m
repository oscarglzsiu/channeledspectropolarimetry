function [intensity] = intensitySCS(stokesincident,phi1,phi2)
% González-Siu L.O., 25/10/2019
% Stokes polarimetry by channeled polarimetry - Simulation

% This function returns the intensity spectrum given the incident SOP and
% the retardances of the SCS (phi1, phi2).
intensity = 0.5*stokesincident(1,:) ...
    + 0.5*stokesincident(2,:).*cos(phi2) ...
    + 0.5*stokesincident(3,:).*sin(phi2).*sin(phi1) ...
    - 0.5*stokesincident(4,:).*sin(phi2).*cos(phi1);
end