# channeledspectropolarimetry
Matlab / Octave library for simulation of passive and active polarimeters with spectral channeling, and extraction of the Stokes vector and the Mueller matrix

This repository contains the matlab functions I wrote during my Master's research on Stokes and Mueller matrix polarimeters with spectral channeling (SCS, and MMCS, respectively). The polarization state analizer (PSA) is composed of two thick birefringent retarders (R1, and R2) followed by a horizontal linear polarizer and a difractive spectrometer. For the MMCS, the PSG is a mirror image of the PSA.

For the SCS, three (3) configurations are considered, based on the orientation of the retarders:
* theta(R1) = 0°, theta(R2) = 45°
* theta(R1) = 0°, theta(R2) = -45°
* theta(R1) = 90°, theta(R2) = -45°

An example (scs_example.m) is provided to explain the extraction process of the Stokes vector as a function of wavelength.

Project publications

[1] González Siu, Luis Oscar, sustentante   Analysis of channeled stokes and mueller matrix polarimeters /   2021

[2] Luis Oscar González-Siu and Neil C. Bruce, “Error analysis of channeled Stokes polarimeters.” Appl. Opt. 60, 4511-4518 (2021) https://www.osapublishing.org/ao/abstract.cfm?uri=ao-60-16-4511

[3] Luis Oscar González-Siu and Neil C. Bruce, “Analysis of experimental errors in Mueller matrix channeled polarimeters.” Appl. Opt. 60, 5456-5464 (2021) https://www.osapublishing.org/ao/abstract.cfm?uri=ao-60-18-5456

Recommended bibliography

[4] Nordsieck, Kenneth H. 1974. “A Simple Polarimetric System for the Lick Observatory Image-Tube Scanner.” Publications of the Astronomical Society of the Pacific 86 (June 1974): 324. https://doi.org/10.1086/129610

[5] Oka, Kazuhiko, and Takayuki Kato. 1999. “Spectroscopic Polarimetry with a Channeled Spectrum.” Optics Letters 24 (21): 1475. https://doi.org/10.1364/OL.24.001475

[6] Iannarilli, Jr., Frank J., Stephen H. Jones, Herman E. Scott, and Paul L. Kebabian. 1999. “Polarimetric-Spectral Intensity Modulation (P-SIM): Enabling Simultaneous Hyperspectral and Polarimetric Imaging.” In Infrared Technology and Applications XXV, edited by Bjorn F Andresen and Marija Strojnik, 3698:474. SPIE. https://doi.org/10.1117/12.354549

[7] Andrey S. Alenin, J. Scott Tyo, “Task-specific snapshot Mueller matrix channeled spectropolarimeter optimization.” Proc. SPIE 8364, Polarization: Measurement, Analysis, and Remote Sensing X, 836402 (8 June 2012); https://doi.org/10.1117/12.921911

[8] Andrey S. Alenin and J. Scott Tyo, “Generalized channeled polarimetry.” J. Opt. Soc. Am. A 31, 1013-1022 (2014) https://doi.org/10.1364/JOSAA.31.001013

