% =========================================================================
%  Project: Cryo-FLM Phase Retrieval and Localization Toolboxes
%
%  Description:
%    This script performs 3D localization for biplane fluorescence datasets 
%    using pre-generated PSF models and CUDA-accelerated MLE fitting.
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


clc;
close all;

%% ------------------------------------------------------------------------
% Add Required Paths
% -------------------------------------------------------------------------
gitpath = '.\Support\Github';
addpath(fullfile(gitpath, 'helpers'));
addpath('.\Support\tag\RC1\PSF Toolbox');

%% ------------------------------------------------------------------------
% Initialize PRPSF Object and Set Parameters
% -------------------------------------------------------------------------
disp('Generating pupil function and estimating aberration using Phase Retrieval');

probj = PRPSF();

% CCD camera properties
probj.CCDoffset = 100;
probj.Gain = 7;

% Optical system configuration
probj.PRstruct.NA = 0.9;
probj.PRstruct.Lambda = 0.515;  % Wavelength in µm
probj.PRstruct.RefractiveIndex = 1;
probj.Pixelsize = 0.130;        % µm/pixel

% Phase retrieval and PSF settings
probj.PSFsize = 256;
probj.SubroiSize = 40;
probj.ZernikeorderN = 7;
probj.PSF_size_forMLE = 16;
probj.resizeFactor = 1;

%% ------------------------------------------------------------------------
% Set Z positions (in microns) for the focal stack
% -------------------------------------------------------------------------
Zpos_plane1 = -1:0.2:1; %for rt
% Zpos_plane1 = -1.4:0.2:1.4; %for cryoT
probj.Zpos = Zpos_plane1;
probj.Zindstart = 1;
probj.Zindend = numel(Zpos_plane1);
probj.Zindstep = 1;

%% ------------------------------------------------------------------------
% Load Bead Data and Noise Data
% -------------------------------------------------------------------------
% Replace ims_Ztrue_plane1_s and ims_Ztrue_plane1_n with your actual variables
input_beadData = ims_Ztrue_plane1_s;   % Signal
input_noiseData = ims_Ztrue_plane1_n;  % Background

probj.BeadData = input_beadData;
probj.NoiseData = input_noiseData;

%% ------------------------------------------------------------------------
% Set Initial Phase Coefficients
% -------------------------------------------------------------------------
num_terms=(probj.ZernikeorderN + 1) * (probj.ZernikeorderN + 2) / 2+2; %aberrations and Intensity&bg
probj.phase_coefficients = zeros(1, num_terms);
probj.phase_coefficients(end-1) = 0.2;
probj.phase_coefficients(end) = 0.1;


%% ------------------------------------------------------------------------
% Preprocessing Steps
% -------------------------------------------------------------------------
probj.prepdata();          % Normalize and format input data
probj.precomputeParam();   % Precompute Fourier and Zernike bases
probj.findXYshift();       % Estimate XY shift between frames
probj.datapreprocess();    % Apply shift correction and cropping

%% ------------------------------------------------------------------------
% Generate Model PSF and Estimate Initial Phase
% -------------------------------------------------------------------------
probj.genMpsf();
probj.estdata();


% Estimate PSF width (sigma) using the updated phase model
probj.precalsigma();
disp(['Calculated Gaussian blur kernel size: [', num2str(probj.gBlur), ']']);

% % Set Gaussian blur parameter manually
% probj.gBlur = 1; %for cryoT beads, use 1.1 for more diffused pattern

%% ------------------------------------------------------------------------
% Output Directory Configuration
% -------------------------------------------------------------------------
save_dir = 'F:\cryoFLM\codes\MLE_publish_version\1. Phase Retrieval\3. Test Data\RT Temperature Bead Stack\';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    disp('Folder created successfully!');
else
    disp('Folder already exists.');
end

%% ------------------------------------------------------------------------
% Phase Retrieval Loop (Iterative Optimization)
% -------------------------------------------------------------------------
max_iterations = 1;

for iteration = 1:max_iterations
    disp(['Iteration ', num2str(iteration), ' of ', num2str(max_iterations)]);

    % Clear previous iteration metrics
    probj.ncc_values = [];
    probj.mse_values = [];
    probj.loglikelihood_values = [];

    % Perform phase retrieval update
    probj.phaseretrieve();

    % Save intermediate results
    save(fullfile(save_dir, ['PR-PSF_iteration_', num2str(iteration), '.mat']), 'probj');

    % Optional: Generate diagnostic plots
    % probj.generate_pm();
    % probj.genPRfigs('PSF');
    % probj.genPRfigs('pupil');
    % probj.plot_similarity(probj.ncc_values, 'Normalized Cross Correlation');
    % probj.plot_similarity(probj.mse_values, 'Mean Square Estimation');
    % probj.plot_similarity(probj.loglikelihood_values, 'Likelihood');
end

disp('Phase-retrieval process completed.');
