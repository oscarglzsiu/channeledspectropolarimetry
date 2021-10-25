% Stokes channeled spectropolarimetry (SCS)
% by Gonzalez-Siu, Luis Oscar (15/09/2021)

% This example considers a Stokes Channeled Spectropolarimeter (SCS)
% consisting of two thick birefringent retarders, R1 and R2, with fast axes
% orientations at 0 and 45 degrees, respectively; next, there is a
% horizontal linear polarizer, with its transmission axis defining the 0
% degrees reference, and a dispersive spectrometer.
% The retarders are made of quartz.

% Clear the workspace memory and close all figures for a clean restart.
% Specially if you want to test different configurations or datasets.
close all;
clear;
clc;

%% 1. Acquire irradiance spectrum I(lambda) ------------------------------
% Including spectrometer details (wavelength bandwidth lambda in nm, No. of
% pixels N_lambda, bandwidth resolution*)
% *resolution of the spectrometer is not necessary to extract the Stokes
% parameters, but it has a high impact on the extraction performance
% Is the light source profile necessary? No, it can be a reference for S0
% extraction, but through normalization, it looses relevance. Although, it
% is used for more realistic simulations.
% The irradiance spectrum must be already recorded. This code will either
% read a data file or simulate it.

% Dim(lambda, and I_lambda) = [1,N_lambda]

% 1.a Read (or simulate) the irradiance spectrum
% 1.b If available, read the extracted Stokes parameters. For simulations,
% calculate the incident Stokes vector.
% 1.c Set SCS parameters
%   dv        - Configuration vector, nominal thickness ratio, called local
%               retardance factors
%   dzero     - Common thickness factor, called global retardance factor
%   d         - Thickness vector. Measured thickness of the retarders
%   crystalid - Retarder crystal ID (crystalid = 1 quartz, 2 calcite)

% SIMULATION ---------------------------------------------------------
% Set the wavelength resolution and shortest wavelength in nm
deltalambda = 0.2;
lambdamin = 543;
lambdamax = 680;
lambda = lambdamin:deltalambda:lambdamax;
N_lambda = length(lambda);
% Light source spectrum
lightsourcefile = ...
  'lightsourceprofile/SLS202L_Spectrum Thorlabs Tungsten.xlsx';
lightsource = readtable(lightsourcefile,'Range','C2:D20087');
tempx = lightsource{:,1}; % wavelength (lambda) of the source in [nm]
tempx = tempx(~isnan(tempx));
tempy = lightsource{:,2}; % intensity (I) of the source in [a.u.]
tempy = tempy(~isnan(tempy));
lightsourceAmp = 10; % arbitrary light source amplitude
% light source intensity interpolated to the spectrometer lambda values
lightsourceintensity = lightsourceAmp*interp1(tempx,tempy,lambda);
% Incident SOP simulation
stokesINid = 3;
StokesIN = FStokes_InitialValues(N_lambda,stokesINid).*...
  lightsourceintensity;
% StokesIN = [1,1,1,1]'.*ones(1,N_lambda); %.*lightsourceintensity;
StokesINnorm = StokesIN(2:end,:)./StokesIN(1,:);
crystalid = [1,1];
% Size of acquired crystals from Newlight Photonics Inc.
% 5, 8, 10, 15, 20 mm
dv = [3 2];
dzero = 1e-3;
d = dv*dzero;
I_lambda = intensitySCS(StokesIN,...
  retardance(birefringenceswitch(crystalid(1),lambda),d(1),lambda),...
  retardance(birefringenceswitch(crystalid(2),lambda),d(2),lambda));

numframes = 10;

%% 2. Extraction of the Stokes parameters --------------------------------
% We have two main branches:
% a. Mueller Stokes formalism
% We use a characteristic matrix of the system Wp, using the first rows of
% the matrices W(lambda)
% b. Channeled polarimetry
% We take advantage of the modulation of the spectrum to analyze the
% channels formed in the Fourier space

% Parameters already set in previous section 1

%% 2.a Extraction of S(lambda) through the characteristic matrix Wp ------
% Calculate the characteristic matrix W(lambda) and cut it into sections,
% hence obtaining a set of sub-matrices Wp representative of their
% corresponding frames.
B1_lambda = birefringenceswitch(crystalid(1),lambda);
B2_lambda = birefringenceswitch(crystalid(2),lambda);
[StokesW,StokesWnorm] = extractSCSmatrixWp(lambda,I_lambda,N_lambda,...
                           numframes,d,B1_lambda,B2_lambda);

%% 2.b Extraction of S(sigma) through channeled methods ------------------
% Channel splitting method (chs), taken from Oka and Kato (1999)
% Analytical channel splitting method (achs), taken from Alenin and Tyo
% (2012,2014)

% It is necessary to convert the (mostly) uniformly spaced spectrum in
% wavelength (lambda) to a uniformly spaced spectrum in wavenumber (sigma).
% Hence, the No. of pixels in the sigma-domain is set by the user.
% Consider the separation between channels in the tau-domain (Fourier
% space) is inversely proportional to N_sigma.
N_sigma = 128;

% Convert the spectrum in wavelength to a uniformly spaced spectrum in
% wavenumber.
[sigma,I_sigma,deltasigma] = ...
  wavelength2wavenumber(lambda,I_lambda,N_sigma);
% Reconstruction artifacts, associated with the Fast Fourier Transform
% (FFT) method applied, appear when extracting the Stokes vector due to a
% lack of periodicity in the irradiance measured, which contradicts the FFT
% supposition of a periodic function. Therefore, apodization has to be
% applied to the raw data. Here, we use a Hann window.
omega = hann(N_sigma).';
% omega = Hannwindow(N_sigma);
I_sigma_omega = I_sigma .* omega;
% Autocorrelation function
[Cfunction,magCfunction] = autocorrelation(I_sigma_omega);
% Filter channels --------------------------------------------------------
% Support variable to convert sigma in [m-1] to lambda in [nm]
lambdafromsigma = 1e9 ./ sigma;
% Birefringence B(sigma) and Retardance phi(sigma)
B1_sigma = birefringenceswitch(crystalid(1), lambdafromsigma);
B2_sigma = birefringenceswitch(crystalid(2), lambdafromsigma);
sigma0index = N_sigma * 0.5 + 1;
sigma0 = sigma(sigma0index);                                 % [m-1]
OPD1_ref = B1_sigma(sigma0index)*d(1)*1e3;                   % [mm]
OPD2_ref = B2_sigma(sigma0index)*d(2)*1e3;                   % [mm]
[Qmatrix,OPDrelpos,A0index] = QmatrixSCS(OPD1_ref,OPD2_ref);

phi1_sigma = retardance( B1_sigma, d(1), lambdafromsigma);
phi2_sigma = retardance( B2_sigma, d(2), lambdafromsigma);
% Index (pixel) position of the channels in the tau-domain
% Instead of the real thickness (d), we use the nominal thickness
% (dv*dzero)

% Channels filtration. The second output are the centered channels.
dB1_sigma = diff(B1_sigma) / deltasigma;
dB2_sigma = diff(B2_sigma) / deltasigma;
L1 = d(1)...
  .*(B1_sigma(sigma0index)+dB1_sigma(sigma0index)*sigma0) * 1e6; % [um]
L2 = d(2)...
  .*(B2_sigma(sigma0index)+dB2_sigma(sigma0index)*sigma0) * 1e6; % [um]
Lvec = [L2,L1+L2,L1-L2];                                         % [um]
deltatau = (1/(sigma(end)-sigma(1)))*1e6;                        % [um]
Lpx = round(Lvec./deltatau);                                     % [px]
indextau = sort(sigma0index + [-Lpx,0,Lpx]);                     % [px]

windowsize = min(abs(Lpx));
% windowfilter = Xrect1d(windowsize * 2,windowsize);
plancktaper_eps = 0.1;
windowfilter = Xplancktaper1d(windowsize,plancktaper_eps);

[AChannels,channelfilterM] = filterChannels(Cfunction,N_sigma,sigma0index,...
  Qmatrix,Lpx,windowfilter,OPDrelpos);

% Channel Splitting method
[Stokeschs,Stokeschsnorm] = extractSCS_channelsplit(...
  AChannels,A0index,OPDrelpos,N_sigma,phi1_sigma,phi2_sigma,dv,omega);

% Analytical Channel Splitting method
[Stokesachs,Stokesachsnorm] = extractSCS_analitycalchannelsplit(...
  AChannels,A0index,Qmatrix,OPDrelpos,phi1_sigma,phi2_sigma,dv,omega);

%% 3. Results ------------------------------------------------------------
% Irradiance spectra: I(lambda) and I(sigma), to verify the conversion.
figure, plot(lambda, I_lambda,1e9./sigma,I_sigma,'+');
xlabel('\lambda \mum');
ylabel('I(\lambda)');
legend({'Recorded spectrum','Interpolated'},...
  'Location','South','NumColumns',2);

% Plot the intensity as I(lambda) and I(sigma), where the second is
% obtained through the intensity equation as a function of the
% retardances and the extracted Stokes parameters.
I_chs = intensitySCS(Stokeschs,phi1_sigma,phi2_sigma);
I_achs = intensitySCS(Stokesachs,phi1_sigma,phi2_sigma);
figure, plot(lambda, I_lambda,1e9./sigma,I_sigma,'--');
hold on;
plot(1e9./sigma, I_chs,'^:', 1e9./sigma, I_achs,'v:');
hold off;
xlabel('\lambda [nm]');
ylim([min(I_lambda), max(I_lambda)]);
legend({'I_{IN}','I_{INT}','I_{CHS}','I_{ACHS}'},...
  'Location','South','NumColumns',4);

%% 3.a Extracted Stokes parameters

% Parameters used for results display
ppframe = floor(N_lambda/numframes);  % pixels per frame

figure;
for i = 1:3
subplot(3,2,i*2);
hold on;
plot(lambda,StokesINnorm(i,:));
plot(1./sigma*1e9,[Stokeschsnorm(i,:);Stokesachsnorm(i,:)],'--');
plot(lambda(ppframe*(0:numframes-1)+floor(ppframe/2)),...
  StokesWnorm(i,:),'ks');
for j = 1:numframes
  plot(lambda(ppframe*(j-1)+[1 ppframe]),...
    StokesWnorm(i,j).*ones(1,2),'k-.');
end
hold off;
xlabel('\lambda [nm]');
ylabel(['S_',num2str(i),'/S_0 [-]']);
xlim([lambda(1), lambda(end)]);
ylim([-1.2 1.2]);
end
legend({'Reference','Channel Split*','A. C. Split*','DRM'},...
    'NumColumns',5,'Location','South');

for i = 1:3
subplot(3,2,i*2-1);
hold on;
plot(lambda,StokesIN(i+1,:));
plot(1./sigma*1e9,[Stokeschs(i+1,:);Stokesachs(i+1,:)],'--');
plot(lambda(ppframe*(0:numframes-1)+floor(ppframe/2)),...
  StokesW(i+1,:),'ks');
for j = 1:numframes
  plot(lambda(ppframe*(j-1)+[1 ppframe]),...
    StokesW(i+1,j).*ones(1,2),'k-.');
end
hold off;
xlabel('\lambda [nm]');
ylabel(['S_',num2str(i),' [a.u.]']);
xlim([lambda(1), lambda(end)]);
%   ylim([-1.2 1.2]);
end
legend({'Reference','Channel Split*','A. C. Split*','DRM'},...
    'NumColumns',5,'Location','South');

%% 3.b Other values

% Autocorrelation function C and filtered channels A_i
plotCChannels(N_sigma,sigma,indextau,magCfunction,AChannels);
