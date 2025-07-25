function phaseretrieve(obj)
    % phaseretrieve - Generate pupil function based on phase retrieval
    % algorithm described in the paper.
    % Edited by Lihongjia, 2024/05/07
    
    % Initialize variables
    z = obj.Zpos;
    k_z = compute_axial_frequency(obj);
    R = obj.PSFsize;
    [coeffVecStart_norm, coeffFactor] = normalize_zernike_coeffs(obj);
   
    disp('Beginning optimization...');

    llfun = @(x) calc_nLogLikelihood(obj, obj.Mpsf_subroi, ...
        calc_img_fromCoeffVec(x ./ coeffFactor, z, ...
        obj.ZernikeorderN, k_z, obj.sVec, obj.bgVec, ...
        obj.gBlur, obj.resizeFactor, R, obj.Z.ZM, ...
        obj.PSF_size_forMLE, obj.NA_constrain));

    [lowerLimVec, upperLimVec] = set_constraints(obj);

    options_fine = optimoptions('fmincon', 'MaxFunEvals', 5000, 'TolFun', 1e-16, ...
                            'TolCon', 1e-8, 'InitBarrierParam', 1, ...
                           'DiffMinChange', 1e-8, 'Display', 'iter', ...
                           'Algorithm', 'interior-point');

    [resVec_fine, objVal_fine] = fmincon(llfun, coeffVecStart_norm, [], [], [], [], ...
                                        lowerLimVec, upperLimVec, [], options_fine);


    resVec=resVec_fine;
    objVal=objVal_fine;

    % Check for optimization issues
    check_optimization(resVec, lowerLimVec, upperLimVec, objVal);
    
    % Finalize coefficients and calculate images
    coeffVec_final = resVec ./ coeffFactor;
    [imgsOut_raw, pupil_phase, pupil_mag] = generate_psf(obj,obj.PSF_size_forMLE, coeffVec_final, z, k_z, R);
    draw_and_save_results(obj, imgsOut_raw, obj.Mpsf_subroi);

    coeffVec_final'
    
    % Update object properties
    update_obj_properties(obj, coeffVec_final, pupil_phase, pupil_mag, imgsOut_raw);
end

% Custom output function to log data for each iteration
function stop = outfun(x, optimValues, state)
    stop = false;
    switch state
        case 'iter'
            fprintf('Iteration %d: Objective Value = %f\n', ...
                optimValues.iteration, optimValues.fval);
    end
end

function k_z = compute_axial_frequency(obj)
    % Compute axial frequency
    k_z = sqrt((obj.PRstruct.RefractiveIndex / obj.PRstruct.Lambda)^2 - obj.k_r.^2) .* obj.NA_constrain;
end

function [coeffVecStart_norm, coeffFactor] = normalize_zernike_coeffs(obj)
    % Normalize Zernike coefficients
    nZernikeCoeffs = (obj.ZernikeorderN + 1) * (obj.ZernikeorderN + 2) / 2 -4; 
    obj.coeffVecStart = [0.2*rand(nZernikeCoeffs,1) - 0.1;0.2;0.1];
    coeffFactor = [ones(nZernikeCoeffs, 1) * 100; 10;10];
    coeffVecStart_norm = obj.coeffVecStart .* coeffFactor;
end

function [lowerLimVec, upperLimVec] = set_constraints(obj)
    lowerLimVec = [-10000000 * ones(length(obj.coeffVecStart) - 2, 1); 0;0];
    upperLimVec = [1000000 * ones(length(obj.coeffVecStart) - 2, 1); 5000;5000];
end

function check_optimization(resVec, lowerLimVec, upperLimVec, objVal)
    % Check for optimization issues
    if isempty(objVal)
        disp('WARNING: objVal is empty');
        objVal = 0;
    end
    if any(abs(resVec - lowerLimVec) < 0.1) || any(abs(resVec - upperLimVec) < 0.1)
        disp('WARNING: Optimization hit a boundary.');
    end
end

function [imgsOut_raw, pupil_phase, pupil_mag] = generate_psf(obj,R_measuredpsf, coeffVec_final, z, k_z, R)
    % Generate PSF images based on final coefficients
    [imgsOut_raw, aberrPhase, pupil_mag, pupil_phase] = calc_img_fromCoeffVec(coeffVec_final, z, ...
        obj.ZernikeorderN, k_z, obj.sVec, obj.bgVec, obj.gBlur, obj.resizeFactor, ...
        R, obj.Z.ZM, R_measuredpsf, obj.NA_constrain);
end

function draw_and_save_results(obj, crop_fittedpsf_raw, crop_measuredpsf)
    draw_comparison(crop_fittedpsf_raw, 'fitted', crop_measuredpsf, 'measured');
    % save2tiff(crop_fittedpsf_raw, 'final_fitted_psf_raw.tif');
    % save2tiff(crop_measuredpsf, 'final_measured_psf.tif');
    save('no_aberration_psf.mat', 'crop_fittedpsf_raw');
    save('measured_psf.mat', 'crop_measuredpsf');
end

function update_obj_properties(obj, coeffVec_final, pupil_phase, pupil_mag, imgsOut_raw)
    coeffVec_final_coeff_phase = [0,0,0,0,coeffVec_final'];
    coeffVec_final_forpupil = [0, 0,0,0,coeffVec_final'];
    obj.phase_coefficients = coeffVec_final_coeff_phase;
    obj.phase_coefficients_pupil = coeffVec_final_forpupil;
    obj.PRstruct.Pupil.phase = pupil_phase;
    obj.PRstruct.Pupil.mag = pupil_mag;
    obj.PSFstruct.PRpsf = imgsOut_raw;
end

%function [imgsOut, aberrPhase, pupil_mag, pupil_phase] = calc_img_fromCoeffVec(coeffVec, zVec, ZernikeorderN, k_z, sVec, bgVec, gBlur, resizeFactor, R, ZM, PSF_size_forMLE, NA_constrain)
function [imgsOut, aberrPhase, pupil_mag, pupil_phase] = calc_img_fromCoeffVec(coeffVec, zVec, ZernikeorderN, k_z, sVec, bgVec, gBlur, resizeFactor, R, ZM, PSF_size_forMLE, NA_constrain)
    % Initialize variables
    N = length(zVec);
    Rpupil_mag = zeros(R, R, N); % retrieved pupil function magnitude
    Rpupil_phase = zeros(R, R, N); % retrieved pupil function phase
    pupil_mag = zeros(R, R);

    nZernikeCoeffs = (ZernikeorderN + 1) * (ZernikeorderN + 2) / 2; 
   
    coeffVec_m = zeros(nZernikeCoeffs,1);
    coeffVec_m(1)=1;
    
    % Generate the pupil magnitude using Zernike polynomials
    for k = 1 : nZernikeCoeffs 
        pupil_mag = pupil_mag + ZM(:, :, k) .* coeffVec_m(k);
    end

    coeffVec=coeffVec';
    coeffVec_p=[0, 0,0,0,coeffVec];
    aberrPhase = zeros(R, R);
    for k = 1 : nZernikeCoeffs
        aberrPhase = aberrPhase + ZM(:, :, k) .* coeffVec_p(k);
    end
    aberrPhase = exp(1i .* aberrPhase);
    
    
    % Initialize the output image stack
    imgsOut = nan([PSF_size_forMLE PSF_size_forMLE N]);
    %imgsOut = nan([R R N]);
    for j = 1:N
        % Compute defocus phase and pupil function
        defocus_phase = 2 * pi * 1i * zVec(j) * k_z;   
        pupil_complex = pupil_mag .* exp(defocus_phase) .* aberrPhase;
        
        Rpupil_mag(:, :, j) = abs(pupil_complex);
        Rpupil_phase(:, :, j) = pupil_complex ./ Rpupil_mag(:, :, j) .* exp(-defocus_phase .* 1i);

        % Generate the PSF from the pupil function
        Eimg_OS = abs(fftshift(fft2(pupil_complex)));
        Iimg_full = Eimg_OS.^2;
        Iimg = Iimg_full ./ sum(Iimg_full(:));

        %for loglikelihood
        % Apply signal scaling and background
        if any(sVec)
            Iimg = Iimg * coeffVec(end-1) * sVec(j)*10 ;
            Iimg = Iimg + max(0, coeffVec(end) * bgVec(j) * 10);
            
        end

        % Apply Gaussian blur if needed
        if gBlur > 0
            filterSize = round(6*gBlur + 1);
            h = fspecial('gaussian', [filterSize filterSize], gBlur);
            Iimg = imfilter(Iimg, h);
        end
        
        % Crop the PSF image to the desired size
        Iimg = Iimg(R/2 - PSF_size_forMLE/2 + 1:R/2 + PSF_size_forMLE/2, ...
                    R/2 - PSF_size_forMLE/2 + 1:R/2 + PSF_size_forMLE/2);
        imgsOut(:, :, j) = Iimg;
    end

    % Normalize the pupil phase and magnitude
    pupil_phase = mean(Rpupil_phase, 3);
    pupil_phase = pupil_phase ./ abs(pupil_phase);
    
    Fig3 = mean(Rpupil_mag, 3) .* NA_constrain;
    Fig4 = abs(pupil_phase) .* Fig3;
    Fig4 = Fig4.^2;
    Fig4 = Fig4 ./ sum(Fig4(:));
    pupil_mag = sqrt(Fig4);
end

function draw_comparison(data1, str1, data2, str2)
    % Create a full-screen figure
    h = figure('Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    
    % Number of frames in each dataset
    nFrames_1 = size(data1, 3);
    nFrames_2 = size(data2, 3);

    % Set number of images per row
    imagesPerRow = 13;
    
    % Calculate the number of rows needed
    nRows_1 = ceil(nFrames_1 / imagesPerRow);
    nRows_2 = ceil(nFrames_2 / imagesPerRow);
    
    % Total rows needed for both datasets
    totalRows = nRows_1 + nRows_2;

    % Loop through and display the first set of images
    for i = 1:nFrames_1
        row = floor((i-1) / imagesPerRow) + 1;
        col = mod(i-1, imagesPerRow) + 1;
        subplot(totalRows, imagesPerRow, (row-1) * imagesPerRow + col);
        imshow(data1(:,:,i), []);  % Display each image
        title(sprintf('%s %d', str1, i));  % Add a title to each subplot
    end

    % Loop through and display the second set of images
    for i = 1:nFrames_2
        row = nRows_1 + floor((i-1) / imagesPerRow) + 1;
        col = mod(i-1, imagesPerRow) + 1;
        subplot(totalRows, imagesPerRow, (row-1) * imagesPerRow + col);
        imshow(data2(:,:,i), []);  % Display each image
        title(sprintf('%s %d', str2, i));  % Add a title to each subplot
    end
    
    % Save the figure as a .jpg file
    saveas(h, 'compare.jpg', 'jpg');
    pause(2);
    % Close the figure
    close(h);
end


function save2tiff(imageArray, filename)
    % Convert the image array to uint16
    minVal = min(imageArray(:));
    maxVal = max(imageArray(:));
    imageArray = uint16((imageArray - minVal) / (maxVal - minVal) * 65535);
    
    % Create a Tiff object
    t = Tiff(filename, 'w');
    
    % Loop through each frame and write to the TIFF file
    for k = 1:size(imageArray, 3)
        % Set the necessary TIFF tags
        t.setTag('ImageLength', size(imageArray, 1));
        t.setTag('ImageWidth', size(imageArray, 2));
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('BitsPerSample', 16);
        t.setTag('SamplesPerPixel', 1);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.setTag('Compression', Tiff.Compression.None);
    
        % Write the image data for the current frame
        t.write(imageArray(:, :, k));
    
        % Write a new directory if there are more frames
        if k < size(imageArray, 3)
            t.writeDirectory();
        end
    end
    
    % Close the Tiff object
    t.close();
end

