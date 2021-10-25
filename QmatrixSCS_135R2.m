function [Qmat,OPDrelpos,A0index] = QmatrixSCS_135R2(B1,B2)
% by Gonzalez-Siu, Luis Oscar (15/09/2021)
% This function returns:
% (1) the Q-matrix of an SCS given it's thickness ratio (dv),
% (2) a vector with the relative position of the (positive) channels with
% respect to the DC channel A_0(tau), and
% (3) the index position of A_0(tau).

% For the relativeposition vector consider the retardances:
% [0, phi2, phi1+phi2, phi1-phi2]

% Support variables related to the relative position of the carrier
% frequencies (L_i, \tau_i, OPD) based on birefringence
sptBvec = [B2, B1+B2, B1-B2];
sptBvecmag = abs(sptBvec);
sptBvecsign = sign(sptBvec);

% Support counters for the No. of channels
% Count of channels, including DC channel (+1), excluding conjugate
% channels (A*_i)
% DC channel A_0(tau) index
A0index = length(unique(sptBvecmag)) + 1;    % sptC = 4, or 3 if dv(2,1)
% Total count of channels, including DC channel (A_0) and conjugate
% channels (A*_i)
% sptNC = sptC * 2 - 1;                     % sptNC = 7, or 5 if dv(2,1)

% We calculate the relative position of the tau_i channel (OPD) comparing
% the magnitude of the pseudo-OPD positions (sptBvecmag) with the unique
% values of the same vector (i.e., unique(sptBvecmag)) and then multiplying
% by the corresposding sign.
if A0index ~= 3
OPDrelpos = (find(sptBvecmag == unique(sptBvecmag).') - 3*(0:2)')...
  .* sptBvecsign';
else
OPDrelpos = [1;2;1];
end

% The Q-matrix template is based on that presented by Gonzalez-Siu and
% Bruce (2021).
QmatrixT_template = [...
        0,      0,    0, 0.5,    0,       0,       0;...
        0,      0, 0.25,   0, 0.25,       0,       0;...
   -0.125, +0.125,    0,   0,    0,  +0.125,  -0.125;...
  +0.125i,-0.125i,    0,   0,    0, +0.125i, -0.125i];
% +(L1-L2) (L1+L2)   +L2   L0   -L2 -(L1+L2) -(L1-L2)
%       3       2     1    0    -1       -2       -3

% Q-matrix transpose definition
QmatrixT = zeros(4,2*A0index-1);
QmatrixT(1,A0index) = 0.5;
for i = 1:3
  QmatrixT(:,A0index + OPDrelpos(i)) = ...
    QmatrixT(:,A0index + OPDrelpos(i)) + ...
    QmatrixT_template(:,4 + i);
  QmatrixT(:,A0index - OPDrelpos(i)) = ...
    QmatrixT(:,A0index - OPDrelpos(i)) + ...
    QmatrixT_template(:,4 - i);
end
Qmat = QmatrixT.';

end