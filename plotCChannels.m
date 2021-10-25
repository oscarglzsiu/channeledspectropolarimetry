function [] = plotCChannels(numpixels,sigma,itau,magCfunction,AChannels)
% by Gonzalez-Siu, Luis Oscar (12/08/2021)
% This function plots the autocorrelation function (magnitude) |C(tau)| and
% the filtered channels (magnitude) |A(tau)|.
tau = taudomainpixelsize(numpixels,sigma);
figure
plot(tau, magCfunction);
hold on;
plot(tau,abs(AChannels),'--');
for i = 1:length(itau)
  xline(tau(itau(i)),'--');
end
hold off;
xlabel('\tau \mum');
ylabel('|C(\tau)|');
end