% This is an example script for using the OTFrescale class

% use either PSF_pupil_example.m or PSF_zernike_example.m to generate a
% stack of PSF images
%% create object and set input properties of OTFrescale class
otfobj=OTFrescale();
otfobj.SigmaX=2;
otfobj.SigmaY=2;
otfobj.Pixelsize=0.1;% micron
otfobj.PSFs=psfobj.PSFs;
%% generate OTF rescaled PSF
otfobj.scaleRspace();
% show PSFs
dipshow(otfobj.Modpsfs); % OTF rescaled PSFs
dipshow(psfobj.PSFs);% input PSFs