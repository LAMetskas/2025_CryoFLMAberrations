clc; close all;

%% Add necessary paths
gitpath = '.\Support\Github';
addpath(fullfile(gitpath, 'Live3DSMSN', 'SRsCMOS'));
addpath(fullfile(gitpath, 'CUDAmex'));
addpath(fullfile(gitpath, 'mex'));
addpath(fullfile(gitpath, 'helpers'));

addpath('.\Support\tag\RC1\PSF Toolbox');

%% Initialize PRPSF object
disp('Visualize the pupil and metrics from phase retrieval...');
probj = PRPSF();

% Define imaging system and acquisition parameters
probj.PRstruct.NA           = 0.9;
probj.PRstruct.Lambda       = 0.575; 
probj.PRstruct.RefractiveIndex = 1.0; 
probj.Pixelsize             = 0.130;   % in microns

% Precompute FFTs and constants needed for optimization
probj.precomputeParam();

%% Phase retrieval iterative loop
save_dir = '..\3. Test Data\Cryo Temperature Bead Stack\';

num_iterations=1;
for iteration = 1:num_iterations
    disp(['Iteration ', num2str(iteration), ' of ', num2str(num_iterations)]);

    % Load intermediate results for the current iteration
    data_path = fullfile(save_dir, ['PR-PSF_iteration_', num2str(iteration), '.mat']);
    load(data_path);

    % Extract data
    raw_psf       = probj.Mpsf_subroi;
    fitted_psf    = probj.PSFstruct.PRpsf;
    coefficients  = probj.phase_coefficients;
    pupil_phase   = probj.PRstruct.Pupil.phase;

    % Draw and save pupil phase
    probj.draw_pupil(save_dir, iteration);

    % Save raw vs fitted PSF image (XY projections)
    raw_psf_img_path    = fullfile(save_dir, ['raw_psf_iteration_', num2str(iteration), '.jpg']);
    fitted_psf_img_path = fullfile(save_dir, ['fitted_psf_iteration_', num2str(iteration), '.jpg']);
    draw_colormap_all_frames_together(raw_psf, fitted_psf, raw_psf_img_path, fitted_psf_img_path);

    % Save XZ plane projections
    raw_psfz_img_path    = fullfile(save_dir, ['raw_psfz_iteration_', num2str(iteration), '.jpg']);
    fitted_psfz_img_path = fullfile(save_dir, ['fitted_psfz_iteration_', num2str(iteration), '.jpg']);
    visualize_zx_plane(raw_psf, fitted_psf, raw_psfz_img_path, fitted_psfz_img_path);

    % Save similarity plots
    probj.plot_similarity(probj.ncc_values, fullfile(save_dir, ['NCC Iteration ', num2str(iteration)]),'NCC');
    probj.plot_similarity(probj.mse_values, fullfile(save_dir, ['MSE Iteration ', num2str(iteration)]),'MSE');
    probj.plot_similarity(probj.loglikelihood_values, fullfile(save_dir, ['LogLikelihood Iteration ', num2str(iteration)]),'LogLikelihood');
end

function visualize_zx_plane(stack1, stack2, save_name1, save_name2)
    % Normalize and clip both stacks for visualization
    stack1 = mat2gray(stack1) / max(stack1(:));
    stack2 = mat2gray(stack2) / max(stack2(:));

    % Extract center y-slice
    y_slice = round(size(stack1, 2) / 2);
    xz1 = squeeze(stack1(:, y_slice, :))';
    xz2 = squeeze(stack2(:, y_slice, :))';

    % Clip intensities
    clip = @(img) mat2gray(min(max(img, prctile(img(:), 1)), prctile(img(:), 100)));
    xz1 = clip(xz1);
    xz2 = clip(xz2);

    % Convert to RGB
    rgb1 = ind2rgb(gray2ind(xz1, 256), parula(256));
    rgb2 = ind2rgb(gray2ind(xz2, 256), parula(256));

    % Save images
    save_rgb_image(rgb1, save_name1, 'X-Z Plane-Measured');
    save_rgb_image(rgb2, save_name2, 'X-Z Plane-Fitted');
end

function draw_colormap_all_frames_together(stack1, stack2, save_name1, save_name2)
    % Normalize and clip each frame
    normalize_clip = @(slice) mat2gray(min(max(slice, prctile(slice(:), 1)), prctile(slice(:), 100)));
    
    % Convert stacks to RGB montage
    to_rgb_stack = @(stack) arrayfun(@(i) ...
        ind2rgb(gray2ind(normalize_clip(stack(:,:,i)), 256), parula(256)), ...
        1:size(stack, 3), 'UniformOutput', false);
    
    rgb_cells1 = to_rgb_stack(stack1);
    rgb1 = cat(4, rgb_cells1{:});
    
    rgb_cells2 = to_rgb_stack(stack2);
    rgb2 = cat(4, rgb_cells2{:});

    % Save montages
    save_montage(rgb1, save_name1, 'Measured');
    save_montage(rgb2, save_name2, 'Fitted');
end

function save_rgb_image(rgb_image, filename, title_str)
    fig = figure('Visible', 'off');
    imshow(rgb_image, 'InitialMagnification', 'fit');
    % title(title_str); 
    xlabel('X'); ylabel('Z'); colorbar;
    saveas(fig, filename);
    close(fig);
end

function save_montage(rgb_stack, filename, title_str)
    fig = figure('Visible', 'off');
    montage(rgb_stack, 'Size', [1 NaN]);
    title(title_str);
    saveas(fig, filename);
    close(fig);
end
