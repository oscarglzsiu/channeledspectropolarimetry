function [S] = FStokes_InitialValues(N,op)
% González-Siu L.O., 13/03/2020
% Stokes polarimetry by channeled polarimetry - Simulation

% This function returns a polarization state given a case selection.

psi = 30/180*pi;                                        % [rad]
Sbase = [...
    1,0,0,0;...   1 non-polarized
    1,1,0,0;...   2 horizontal linearly polarized light
    1,0,1,0;...   3 +45 linearly polarized light
    1,0,0,1;...   4 right linearly polarized light
    1,0,1/sqrt(2),1/sqrt(2);...       5
    1,cos(2.*psi),sin(2.*psi),0;...   6
    1,1/sqrt(2),1/sqrt(2),0;...       7
    1,1/sqrt(2),0,1/sqrt(2);...       8
    1,1/sqrt(3),1/sqrt(3),1/sqrt(3);...   9
    1,-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)];  % 10
% S = stokesBase(Sbase(op,:)',N);
S = Sbase(op,:)'.*ones(1,N);
end