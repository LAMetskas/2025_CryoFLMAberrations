% =========================================================================
%  Project: Cryo-FLM Phase Retrieval and Localization Toolboxes
%
%  Description:
%    This script performs biplane fluorescence datasets simulation 
%
%  Copyright (c) 2025 The Huang Lab
%  All Rights Reserved.
%
%  Laboratory:
%    The Huang Lab
%    Weldon School of Biomedical Engineering
%    Purdue University
%    West Lafayette, Indiana, USA
%
%  Author: Hongjia Li
%  Version: 1.0
%  Date: 2025-07-24
% =========================================================================

clc; close all;

%% === Add Required Toolbox ===
addpath('./toolbox/tag/RC1/PSF Toolbox/');

%% === Initial Parameters ===
imsz = 32;                          % Image patch size (pixels)
Numimage = 15;                      % Number of Z-slices
frames_per_z = 100;                 % Frames per Z position
photon_numbers = [2000, 5000, 10000, 20000, 50000, 100000, 200000, 400000]; % Photon counts
bg_num = 400;                       % Background photons (total)
dist_Bi = 1;                        % Biplane separation (µm)
start_z = -0.7;                     % Starting axial position (µm)
step_size = 0.1;                    % Z-step size (µm)
pixel_size = 0.13;                  % Pixel size (µm)

%% === Load Aberration Model ===
% Choose the aberration model from a phase retrieval output.
aberrs = load('../1. Phase Retrieval/3. Test Data/1. Bead Stack RT/Results/PR-PSF_iteration_1.mat');

%% === Main Loop Over Photon Numbers ===
for photon = photon_numbers
    
    %% Create Output Folder
    resultfolder = sprintf("./test_data_simulation/biplanedis%.2f_bg_%d_intensity_%d", ...
        dist_Bi, bg_num, photon);
    if ~exist(resultfolder, 'dir')
        mkdir(resultfolder);
        disp(['Created folder: ', resultfolder]);
    else
        disp(['Folder already exists: ', resultfolder]);
    end

    %% Initialize Image Stacks
    imsp1 = zeros(imsz, imsz, Numimage, frames_per_z);  % Biplane channel 1
    imsp2 = zeros(imsz, imsz, Numimage, frames_per_z);  % Biplane channel 2

    photon_num = photon / 2;  % Half photons per plane
    bg = bg_num / 2;          % Half background per plane

    %% Loop Over Z Positions
    for num = 1:Numimage
        Zpos = start_z + step_size * (num - 1);  % Axial position

        % --- PSF Configuration (Common) ---
        PRstruct = aberrs.probj.PRstruct;
        PRstruct.Lambda = 0.515;                                     % Wavelength (µm)
        PRstruct.Pupil.phase = zeros(128);                           % Initialize pupil phase
        PRstruct.Pupil.mag = zeros(128);                             % Initialize pupil magnitude
        PRstruct.Zernike_phase = aberrs.probj.phase_coefficients(1:end-2); % Zernike phase terms
        PRstruct.Zernike_mag = [1, zeros(1, 35)];                    % Zernike magnitude terms
        PRstruct.SigmaX = 2;                                         % Gaussian sigma X
        PRstruct.SigmaY = 2;                                         % Gaussian sigma Y

        % --- Generate PSF for Plane 1 (Z - dist/2) ---
        psfobj1 = PSF_zernike(PRstruct);
        psfobj1.Xpos = 0; psfobj1.Ypos = 0;
        psfobj1.Zpos = Zpos - dist_Bi / 2;
        psfobj1.Boxsize = imsz;
        psfobj1.Pixelsize = pixel_size;
        psfobj1.PSFsize = 128;
        psfobj1.nMed = 1.33;
        psfobj1.precomputeParam();
        psfobj1.genZernike();
        psfobj1.genPupil();
        psfobj1.genPSF_2();
        psfobj1.scalePSF('normal');

        % --- Generate PSF for Plane 2 (Z + dist/2) ---
        psfobj2 = PSF_zernike(PRstruct);
        psfobj2.Xpos = 0; psfobj2.Ypos = 0;
        psfobj2.Zpos = Zpos + dist_Bi / 2;
        psfobj2.Boxsize = imsz;
        psfobj2.Pixelsize = pixel_size;
        psfobj2.PSFsize = 128;
        psfobj2.nMed = 1.33;
        psfobj2.precomputeParam();
        psfobj2.genZernike();
        psfobj2.genPupil();
        psfobj2.genPSF_2();
        psfobj2.scalePSF('normal');

        %% Generate Frames for Current Z Position
        for frame = 1:frames_per_z
            % Normalize PSFs and scale by photon count
            img1 = (psfobj1.ScaledPSFs / sum(psfobj1.ScaledPSFs(:))) * photon_num;
            img2 = (psfobj2.ScaledPSFs / sum(psfobj2.ScaledPSFs(:))) * photon_num;

            % Add Poisson noise
            noisy_img1 = poissrnd(img1 + bg);
            noisy_img2 = poissrnd(img2 + bg);

            % Store in stacks
            imsp1(:, :, num, frame) = noisy_img1;
            imsp2(:, :, num, frame) = noisy_img2;
        end
    end

    %% Reshape to [x, y, t] where t = Numimage * frames_per_z
    subregion_ch1 = reshape(permute(double(imsp1), [1, 2, 4, 3]), imsz, imsz, []);
    subregion_ch2 = reshape(permute(double(imsp2), [1, 2, 4, 3]), imsz, imsz, []);

    %% Save Dataset
    save(fullfile(resultfolder, 'subregions_beads1.mat'), 'subregion_ch1', 'subregion_ch2');
    disp(['Saved dataset: ', resultfolder]);
end
