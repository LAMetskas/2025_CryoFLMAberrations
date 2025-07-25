% sample code for generating PSF with index mismatch aberraton
% reference : Sheng Liu, Optics Express, 2013

% theta1: angle inside immersion medium
% theta2: angle inside sample medium

nImm = 1.51;               % refractive index of immersion medium
nMed = 1.33;               % refractive index of sample medium
Lambda = 0.6;              % emission wavelength, um 
NA = 1.49;                 % numerical aperture of the objective lens
PSFsize = 128;             % image size of simulated psf
Pixelsize = 0.1;           % pixel size at sample plane, um
stagepos = 2;              % position of the sample stage, um 
z0 = stagepos*nMed/nImm;   % reference z position
zpos = [-z0:0.1:1];        % relative z position of emitter in the sample medium, reference z0, um

%% simulate k space

[X,Y] = meshgrid(-PSFsize/2:PSFsize/2-1,-PSFsize/2:PSFsize/2-1);
Zo = sqrt(X.^2+Y.^2);
scale = PSFsize*Pixelsize;
k_r = Zo./scale;
sin_theta1 = k_r.*Lambda./nImm;
sin_theta2 = nImm./nMed.*sin_theta1;
Cos1 = sqrt(1-sin_theta1.^2);
Cos2 = sqrt(1-sin_theta2.^2);
Freq_max = NA/Lambda;
pupil = k_r<Freq_max;

%% aberration phase from index mismatch
deltaH = z0*nMed.*Cos2 - z0*nImm^2/nMed.*Cos1;
IMMphase = exp(2*pi/Lambda.*deltaH.*1i);

%% generate psf
N = numel(zpos);
psfs = zeros(PSFsize,PSFsize,N);
for ii = 1:N
    defocusMed = exp(2*pi/Lambda*nMed*zpos(ii).*Cos2.*1i);
    pupil_complex = defocusMed.*IMMphase.*pupil;
    psfA = abs(fftshift(fft2(pupil_complex)));
    psfs(:,:,ii) = psfA.^2./R^2;
end
%% IMMphase.*pupil --> the pupil of index mismatch, edited by FX
dipshow(abs(IMMphase.*pupil));
dipshow(angle(IMMphase.*pupil));

psf_xz = permute(squeeze(psfs(PSFsize/2,:,:)),[2,1]);

h = figure('position',[300,300,PSFsize*5,N*5]);
imagesc(psf_xz)
colormap(hot)



