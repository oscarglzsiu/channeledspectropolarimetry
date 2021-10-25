function [Achannels,channelfilterM] = filterChannels...
  (Cfunction,numpixels,sigma0index,Qmatrix,Lpx,windowfilter,OPDrelpos)
% by Gonzalez-Siu, Luis Oscar (15/09/2021)
% This function returns the filtered and centered A channels of an
% autocorrelation function C(tau). This is done using the expected
% positions (indextau) and the Q-matrix specific of the SCS configuration
% (thickness ratio). 

% AChannels [#C x N] where AChannels(i,:) = Channel A_i
% "\bold{C} is a length N vector that lists the composition of each
% channel". (Li2020)

% Initialize channels 'vectors'. Each row is a vector on itself contaning
% a filtered channel.
Achannels = zeros(size(Qmatrix,1),numpixels);
% Filter vector (matrix composed by filters for each channel)
channelfilterM = zeros(size(Achannels));

deltafuncOPD = zeros(1,numpixels);
deltafuncOPD(sigma0index) = 1;
channelfilter = conv(deltafuncOPD,windowfilter,'same');
A0index = ceil(size(Qmatrix,1) * 0.5);
Achannels(A0index,:) = Cfunction .* channelfilter;
channelfilterM(A0index,:) = channelfilter;
for i = 1:length(OPDrelpos)
deltafuncOPD = zeros(1,numpixels);
deltafuncOPD(sigma0index + Lpx(i)) = 1;
% disp(sigma0index + Lpx(i));
channelfilter = conv(deltafuncOPD,windowfilter,'same');
Achannels(A0index + OPDrelpos(i),:) = Cfunction .* channelfilter;
channelfilterM(A0index + OPDrelpos(i),:) = channelfilter;
deltafuncOPD = zeros(1,numpixels);
deltafuncOPD(sigma0index - Lpx(i)) = 1;
% disp(sigma0index - Lpx(i));
channelfilter = conv(deltafuncOPD,windowfilter,'same');
Achannels(A0index - OPDrelpos(i),:) = Cfunction .* channelfilter;
channelfilterM(A0index - OPDrelpos(i),:) = channelfilter;
end

end