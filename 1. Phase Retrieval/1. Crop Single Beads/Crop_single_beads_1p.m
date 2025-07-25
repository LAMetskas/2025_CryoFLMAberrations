% =========================================================================
%  Project: Cryo-FLM Phase Retrieval and Localization Toolboxes
%
%  Description:
%    This script extracts a high-contrast subregion containing a single bead
%    from a 3D (Z-stack) fluorescence image. It is used as the preprocessing
%    step for PSF phase retrieval.
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

% Add DIPimage to MATLAB path (modify the path below to match your installation)
addpath('C:\software\DIPimage 2.9\common\dipimage');
dip_initialise;
dipsetpref('ImageFilePath', 'C:\software\DIPimage 2.9\common\dipimage');
dipsetpref('DefaultMappingMode','lin');
dipsetpref('TrueSize','off');
set(0,'DefaultFigurePaperType','US Letter');
set(0,'DefaultAxesFontName', 'Arial');
set(0,'DefaultAxesFontSize', 15);
clc; close all;

load('..\3. Test Data\1. Bead Stack RT\Beads_stack.mat')

%% Add paths for required functions
gitpath = '..\2. Phase Retrieval\Support\Github';
addpath(fullfile(gitpath, 'helpers'));

%% Load image stack
% Expected shape: [256 x 256  x 41]
% Averaging across channels (if needed): ims_plane1 = squeeze(mean(ims, 3));
ims_plane1 = ims;  % Using full input stack as-is

% Max intensity projection over z-stack (dim 3) to find bright regions
img1_max = max(ims_plane1, [], 3);

%% Select subregion based on brightest bead
boxsize = 64;  % Size of the initial FOV for bead picking
N = 1;         % Number of beads to pick (choose the brightest)

% Pick bead: return FOV and coordinates of brightest bead
[fov, centers] = pickbead(img1_max, boxsize, N, 'stack', []);

% Extract center coordinates
row = centers(2);
col = centers(1);
fprintf('Selected center: row = %d, col = %d\n', row, col);

%% Extract subregion for phase retrieval
sz = 32;  % Half-size of the cropped subregion
% Final size: 64 × 64 × num_frames

% Extract subregion centered on (row, col)
subregion_ch1_s = ims_plane1(row-sz+1:row+sz, col-sz+1:col+sz, :);

% Select all frames (z-slices)
ims_Ztrue_plane1_s = subregion_ch1_s(:, :, 1:end);

% Display with linear scaling
dipshow(ims_Ztrue_plane1_s, 'lin');

% Save subregion
save('subregion_ch1_s.mat', 'subregion_ch1_s');
