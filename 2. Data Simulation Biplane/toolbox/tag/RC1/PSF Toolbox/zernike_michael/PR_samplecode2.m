clear all
close all

% load in the PRtest struct and the PSF Data
load('PRtestSample.mat')
% explanation of PRtest struct:
%
% PRtest.ZernikeOrder : Zernike polynomial index n, defined in Wyant's
%                       paper                       
% PRtest.Lambda: wavelength of fluorescence emission from beads
% PRtest.NA: numerical aperture of the objective length
% PRtest.PixelSizeCCD: pixel size of the CCD chip
% PRtest.Magnify: magification of the microscope
% PRtest.OutputPsfSize: output size of the generated PSF
% PRtest.CN_complex: complex Zernike coefficients of the Zernike polymials
%                    upto nth order
% PRtest.RefractiveIndex: refractive index of the immersion medium of the
%                         objective lens
% PRtest.Zstack: a vector of the z posistions of the measured PSF Data
% Data: measured PSF data, pixel values are converted to the photon counts

n = 1000;
tic
for i = 1 : n
[PSF] = Generate_PSF(PRtest);
end
toc

% shows measured PSF
if n == 1
   dipshow(Data)
   % shows generated PSF, PSF at each z position is normalized, meaning the sum
   % over all pixels is 1.
   dipshow(PSF)
end

e = Err(Data, PSF);
