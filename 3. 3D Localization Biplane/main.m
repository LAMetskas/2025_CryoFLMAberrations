% =========================================================================
%  Project: Cryo-FLM Phase Retrieval and Localization Toolboxes
%
%  Description:
%    This script performs phase retrieval based on MLE
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

%% ========================== 1. PARAMETERS ===============================
% --- Paths ---
addpath('./toolbox/tag/RC1/PSF Toolbox/');

% --- Input / Output ---
data_file = ['..\2. Data Simulation Biplane\test_data_simulation\' ...
    'biplanedis1.00_bg_400_intensity_10000\subregions_beads1.mat'];
aberration_file = ['../1. Phase Retrieval/3. Test Data/1. Bead Stack RT/' ...
    'Results/PR-PSF_iteration_1.mat'];
save_root = fullfile('.', 'loc_results\');

% --- Bead & Model ---
bead_idx = 1;               % Bead index
model_idx = 1;              % Model index

% --- Initial Values ---
I_init = 10000 / 2;        % Initial intensity guess per channel
bg_init = 400 / 2;          % Initial background guess per channel

% --- XYZ Guessing ---
x_guess = 16.5;             % Initial x-position (for imsz=32)
y_guess = 16.5;             % Initial y-position (for imsz=32)
z_values = -0.7:0.1:0.7;    % Z positions for initial guess
frames_per_z = 100;         % Frames per Z position

% --- PSF Parameters ---
lambda_override = 0.515;    % Wavelength override (µm)
psf_size = 128;             % PSF grid size
dist_biplane = 1;           % Biplane separation (µm)
bin = 4;                    % Binning factor
Nzs = 601;                  % Number of Z-slices

% --- CUDA Localization ---
iterateN = 400;             % Number of iterations
lambda_reg = 0;             % Regularization parameter

% --- Output Directory ---
save_dir = fullfile(save_root, sprintf('intensity%d_init', I_init * 2));
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

%% ========================== 2. LOAD DATA ================================
load(data_file, 'subregion_ch1', 'subregion_ch2');  % Load bead data
aberrs = load(aberration_file);                     % Load aberration model

%% ========================== 3. CONSTRUCT PSF OBJECT =====================
oprobj = aberrs.probj;
PRstruct = oprobj.PRstruct;
PRstruct.Lambda = lambda_override;                  % Override wavelength

bxsz = size(subregion_ch1, 1);                      % Subregion size
pxsz = oprobj.Pixelsize;                            % Pixel size

psfobj = PSF_zernike_cryoFM(PRstruct);
psfobj.ZernikeorderN = oprobj.ZernikeorderN;
psfobj.Boxsize = bxsz;
psfobj.Pixelsize = pxsz;
psfobj.PSFsize = psf_size;
psfobj.nMed = PRstruct.RefractiveIndex;
psfobj.PRstruct.SigmaX = 2;
psfobj.PRstruct.SigmaY = 2;

% Remove bias terms from Zernike phase coefficients
psfobj.phase_coefficients = oprobj.phase_coefficients(1:end-2);
psfobj.precomputeParam();
psfobj.genZernike();
psfobj.genPupil();

%% ========================== 4. GENERATE BIPLANE PSFS ====================
[samplepsf, startx, starty, startz, dz, dx] = ...
    gensamplepsf_biplane(psfobj, pxsz, psf_size * bin, bxsz, bin, Nzs, dist_biplane);

% Convert to CUDA-friendly format
samplepsf_cuda = [];
samplepsf_4d = [];
for ii = 1:2
    reshaped = permute(reshape(samplepsf{ii}, bxsz * bin, bxsz * bin, Nzs), [2, 1, 3, 4]);
    samplepsf_cuda = cat(3, samplepsf_cuda, reshaped);
    samplepsf_4d = cat(4, samplepsf_4d, reshaped);
end
samplepsf_cuda = single(samplepsf_cuda);

%% ========================== 5. SPLINE FITTING ===========================
st1 = genpsfstruct(samplepsf_4d(:, :, :, 1), dx, dz, 'matrix');
st2 = genpsfstruct(samplepsf_4d(:, :, :, 2), dx, dz, 'matrix');
st = catstruct(st1, st2, 3);

%% ========================== 6. PREPARE DATA =============================
Nfit = size(subregion_ch1, 3);
data = cat(4, subregion_ch1, subregion_ch2);  % [X Y Frame Channel]
data_cuda = single(data(:));

coords_cuda = single(repmat([1; 1; 0], 1, Nfit));
coords_cuda = coords_cuda(:);

z_guess = repelem(z_values, frames_per_z)';

xtmp = ones(1, Nfit) * x_guess;
ytmp = ones(1, Nfit) * y_guess;

I_next = repmat(I_init, Nfit, 2);   % Initial intensity
bg_next = repmat(bg_init, Nfit, 2); % Initial background

x0 = cat(2, xtmp', ytmp', z_guess, I_next, bg_next);
x0i = single(reshape(x0', Nfit * 7, 1));  % Flattened initial guess

%% ========================== 7. CUDA LOCALIZATION ========================
[P_cuda, ~, crlb, ~, ~] = mex_loc_spline_biplane_32x( ...
    data_cuda, coords_cuda, samplepsf_cuda, ...
    dx, dz, startx, starty, startz, ...
    iterateN, Nfit, lambda_reg, x0i, ...
    single(st.Fx), single(st.Fy), single(st.Fz), ...
    single(st.Fxy), single(st.Fxz), single(st.Fyz), single(st.Fxyz) ...
);
clear mex_loc_spline_biplane_32x;

%% ========================== 8. POST-PROCESSING ==========================
PcudaM = reshape(P_cuda, 7, Nfit)';  % [x y z I1 I2 bg1 bg2]
loc_x = PcudaM(:, 1);
loc_y = PcudaM(:, 2);
loc_z = PcudaM(:, 3);

crlbM = reshape(crlb, 7, Nfit)';
crlb_x = crlbM(:, 1);
crlb_y = crlbM(:, 2);
crlb_z = crlbM(:, 3);

% Save localization results
save(fullfile(save_dir, sprintf('beads%d_model%d_xyzpos.mat', bead_idx, model_idx)), ...
     'loc_x', 'loc_y', 'loc_z');

% Plot z-localization
f_zlocs = figure;
plot(loc_z, '.');
ylabel('Z Localization'); xlabel('Frame');
ylim([-0.8, 0.8]); yticks(-0.8:0.2:0.8);

saveas(f_zlocs, fullfile(save_dir, sprintf('beads%d_model%d_zpos.jpg', bead_idx, model_idx)));
saveas(f_zlocs, fullfile(save_dir, sprintf('beads%d_model%d_zpos.fig', bead_idx, model_idx)));

disp("Finished!");
