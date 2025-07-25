% =========================================================================
%  Project: Cryo-FLM Phase Retrieval and Localization Toolboxes
%
%  Description:
%    This script extracts a high-contrast subregion containing a background
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

%% Add necessary paths
% Update this pathto match your local git repository
gitpath = '..\2. Phase Retrieval\Support\Github';
addpath(fullfile(gitpath, 'helpers'));

%% Load or define the 4D input image stack
% Input format expected: [256 × 256 × n]
% If already averaged or in correct shape, directly use:
ims_plane1 = ims;

% Generate max projection along Z for localization
img1_max = max(ims_plane1, [], 3);

%% Subregion selection
% Parameters
boxsize = 64;  % Size of the initial crop around bead
N = 1;         % Number of beads to pick

% Select the brightest bead region from the image
% Output: fov - cropped image, centers - coordinates of bead center
[fov, centers] = pickbead(img1_max, boxsize, N, 'stack', []);

% Extract center coordinates
row = centers(2);
col = centers(1);

%% Extract high-contrast subregion around selected bead
% Define half-size of target subregion
sz = 32; % Final subregion will be 64×64 in xy

% Crop subregion centered at selected coordinates
% Output size: [64 × 64 × num_frames]
subregion_ch1_n = ims_plane1(row-sz+1:row+sz, col-sz+1:col+sz, :);

% Optionally reduce or select frames (here we use all)
ims_Ztrue_plane1_n = subregion_ch1_n(:, :, 1:1:end);

% Visualize result
dipshow(ims_Ztrue_plane1_n, 'lin');

%% Save the subregion
save('subregion_ch1_n.mat', 'subregion_ch1_n');
