function W = mueller_SCSW_135R2(phi1,phi2)
% This approach considers a non-polarized light source.
W = 0.5 * ones(length(phi1),4);
W(:,2:4) = W(:,2:4) .* [cos(phi2);...
                        -sin(phi1).*sin(phi2);...
                        +cos(phi1).*sin(phi2)].';
% Depricated W-matrix calculation
%  Retarders orientation thetaR1 = 0, thetaR2 = pi/4
% MR1(1:4,1:4,1:numpixels) = mueller_retarder(phi1, 0);
% MR2(1:4,1:4,1:numpixels) = mueller_retarder(phi2, -pi/4);
%  Horizontal linear polarizer
% MP = mueller_polarizer(0);
% W = zeros(4,4,numpixels);
% for i=1:numpixels
%     W(:,:,i) = MP*MR2(:,:,i)*MR1(:,:,i);
% end
end