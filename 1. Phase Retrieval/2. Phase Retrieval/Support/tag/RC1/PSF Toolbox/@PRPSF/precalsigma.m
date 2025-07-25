function precalsigma(obj)
% precalsigma_new - Estimate Gaussian blur (gBlur) by optimizing similarity
% between simulated and measured PSFs via normalized cross-correlation (NCC).
% Updated by Hongjia Li, 2024-05-07

    % Extract parameters
    z = obj.Zpos;
    R = obj.PSFsize;
    coeffVec_final = obj.phase_coefficients;

    % Calculate axial spatial frequency component with NA mask
    k_z = sqrt((obj.PRstruct.RefractiveIndex / obj.PRstruct.Lambda)^2 - obj.k_r.^2) ...
          .* obj.NA_constrain;

    % Get the measured PSF (cropped)
    zdim = size(obj.Mpsf_subroi, 3);
    middle= round(zdim/2);
    if middle-2 >= 1 && middle+2 <= zdim
        crop_measuredpsf = obj.Mpsf_subroi(:,:,middle-2:middle+2);
    else
        error('Index out of bounds: Adjust the value of "middle".');
    end

    % Simulate PSF stack from current Zernike coefficients
    [imgsOut_raw, ~, ~] = generate_psf(obj, obj.PSF_size_forMLE, coeffVec_final, z, k_z, R);
    zdim = size(imgsOut_raw, 3);
    middle= round(zdim/2);
    if middle-2 >= 1 && middle+2 <= zdim
        imgsOut_raw = imgsOut_raw(:,:,middle-2:middle+2);
    else
        error('Index out of bounds: Adjust the value of "middle".');
    end

    % Estimate optimal Gaussian blur to match measured PSFs
    obj.gBlur = cal_gaussian_ncc(imgsOut_raw, crop_measuredpsf);
end

function bestSigma = cal_gaussian_ncc(sourceStack, targetStack)
% cal_gaussian_ncc - Find the optimal Gaussian sigma to align simulated
% and measured PSFs using average slice-wise normalized cross-correlation (NCC).

    % sigmaValues = 0.5:0.1:2.0; % Search range for sigma
    sigmaValues = 0.5:0.1:1.5; % Search range for sigma
    nccValues = zeros(size(sigmaValues));  % Preallocate NCC results
    maxNCC = -inf;
    bestSigma = 0;

    for i = 1:length(sigmaValues)
        sigma = sigmaValues(i);
        filterSize = round(6 * sigma + 1);  % Typical Gaussian kernel size

        % Apply 2D Gaussian blur to each slice of simulated PSF stack
        blurredStack = nan(size(sourceStack));
        h = fspecial('gaussian', [filterSize filterSize], sigma);
        for j = 1:size(sourceStack, 3)
            blurredStack(:, :, j) = imfilter(sourceStack(:, :, j), h, 'replicate');
        end

        % Compute NCC per slice
        ncc_slices = zeros(size(sourceStack, 3), 1);
        for k = 1:size(sourceStack, 3)
            ncc_slices(k) = corr2(blurredStack(:, :, k), targetStack(:, :, k));
        end
        nccValue = mean(ncc_slices);  % Average NCC

        nccValues(i) = nccValue;

        % Track best sigma
        if nccValue > maxNCC
            maxNCC = nccValue;
            bestSigma = sigma;
        end
    end

    disp('NCC values for each sigma:');
    disp(nccValues);
end

function [imgsOut_raw, pupil_phase, pupil_mag] = generate_psf(obj, PSF_size_forMLE, coeffVec_final, z, k_z, R)
% generate_psf - Generate simulated PSFs from Zernike coefficients

    [imgsOut_raw, ~, pupil_mag, pupil_phase] = calc_img_fromCoeffVec( ...
        coeffVec_final, z, obj.ZernikeorderN, k_z, obj.sVec, obj.bgVec, ...
        R, obj.Z.ZM, PSF_size_forMLE, obj.NA_constrain ...
    );
end

function [imgsOut, aberrPhase, pupil_mag, pupil_phase] = calc_img_fromCoeffVec( ...
    coeffVec, zVec, nZernikeCoeffs, k_z, sVec, bgVec, R, ZM, PSF_size_forMLE, NA_constrain)
% calc_img_fromCoeffVec - Compute simulated PSF stack from Zernike coefficients

    N = length(zVec);
    imgsOut = nan([PSF_size_forMLE PSF_size_forMLE N]);
    Rpupil_mag = zeros(R, R, N);
    Rpupil_phase = zeros(R, R, N);
    pupil_mag = zeros(R, R);

    % Amplitude part: magnitude = linear combination of Zernike modes
    coeffVec_m = zeros(nZernikeCoeffs+1, 1);
    coeffVec_m(1) = 1;
    for k = 1:nZernikeCoeffs+1
        pupil_mag = pupil_mag + ZM(:, :, k) .* coeffVec_m(k);
    end

    % Phase part: complex exponential of Zernike phase
    coeffVec_p = [0; coeffVec(:)];
    aberrPhase = zeros(R, R);
    for k = 1:nZernikeCoeffs+1
        aberrPhase = aberrPhase + ZM(:, :, k) .* coeffVec_p(k);
    end
    aberrPhase = exp(1i .* aberrPhase);

    % Generate PSF at each z-position
    for j = 1:N
        defocus_phase = 2 * pi * 1i * zVec(j) * k_z;
        pupil_complex = pupil_mag .* exp(defocus_phase) .* aberrPhase;

        % Store magnitude and normalized phase
        Rpupil_mag(:, :, j) = abs(pupil_complex);
        Rpupil_phase(:, :, j) = pupil_complex ./ Rpupil_mag(:, :, j) .* exp(-1i * defocus_phase);

        % Generate image by taking squared Fourier magnitude
        Eimg = abs(fftshift(fft2(pupil_complex)));
        Iimg = Eimg.^2;
        Iimg = Iimg ./ sum(Iimg(:));

        % Apply signal scaling and background
        if any(sVec)
            Iimg = Iimg * coeffVec(end-1) * sVec(j) * 10;
            Iimg = Iimg + max(0, coeffVec(end) * bgVec(j) * 10);
        end

        % Crop to desired size
        Iimg = Iimg(R/2 - PSF_size_forMLE/2 + 1:R/2 + PSF_size_forMLE/2, ...
                    R/2 - PSF_size_forMLE/2 + 1:R/2 + PSF_size_forMLE/2);
        imgsOut(:, :, j) = Iimg;
    end

    % Normalize pupil phase across Z
    pupil_phase = mean(Rpupil_phase, 3);
    pupil_phase = pupil_phase ./ abs(pupil_phase);

    % Normalize pupil magnitude
    Fig3 = mean(Rpupil_mag, 3) .* NA_constrain;
    Fig4 = abs(pupil_phase) .* Fig3;
    Fig4 = Fig4.^2;
    Fig4 = Fig4 ./ sum(Fig4(:));
    pupil_mag = sqrt(Fig4);
end
