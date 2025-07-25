function draw_pupil(obj, save_dir, iteration)
    % DRAW_PUPIL - Generate and visualize only RMS-normalized phase pupil function
    %
    % Inputs:
    %   obj        - Object containing Zernike polynomials and coefficients
    %   save_dir   - Directory where output images are saved
    %   iteration  - Iteration index for saving figure filenames
    %
    % This function:
    %   1. Computes the RMS-normalized phase pupil function using Zernike polynomials.
    %   2. Visualizes only the RMS-normalized phase.
    %   3. Saves the figure as an image.

    %% === Parameters and Initialization ===
    R = 256;  % Resolution of the pupil grid
    n = ceil((obj.ZernikeorderN + 1) * (obj.ZernikeorderN + 2) / 2);  
    % Number of Zernike coefficients based on order N

    pupil_phasenorm = zeros(R, R);  % Initialize phase pupil matrix

    %% === Compute Phase Pupil Function ===
    for k = 1:n
        pupil_phasenorm = pupil_phasenorm + obj.Z.ZM(:, :, k) .* obj.phase_coefficients(k);
    end

    %% === RMS Phase Normalization ===
    valid_pupil = obj.Z.ZM(:, :, 1) > 0;  % First Zernike mode defines aperture
    phase_values = pupil_phasenorm(valid_pupil);
    
    mean_phase = mean(phase_values);
    rms_phase = sqrt(mean((phase_values - mean_phase).^2));
    pupil_phase_rms = (pupil_phasenorm - mean_phase) / rms_phase;

    %% === Store Result ===
    obj.PRstruct.Fittedpupil.rms_phase = pupil_phase_rms;

    %% === Visualization (Only RMS Phase) ===
    fig = figure('Color', [1, 1, 1], ...
                 'Name', 'RMS-Normalized Pupil Phase', ...
                 'Units', 'normalized', ...
                 'Position', [0.3, 0.3, 0.5, 0.5]);

    RC = 128;   % Center of pupil
    Rsub = 63;  % Subregion size

    imagesc(double(pupil_phase_rms(RC - Rsub : RC + Rsub, RC - Rsub : RC + Rsub)));
    colormap(parula);
    colorbar;
    clim([-5, 5]);  % Adjust phase range as needed
    axis equal off;
    title('RMS-Normalized Phase', 'FontSize', 20, 'FontWeight', 'bold');

    %% === Save Figure ===
    saveas(fig, fullfile(save_dir, ['pupil_rms_phase_iteration_', num2str(iteration), '.jpg']));
    % Figure remains open (not closed)
end
