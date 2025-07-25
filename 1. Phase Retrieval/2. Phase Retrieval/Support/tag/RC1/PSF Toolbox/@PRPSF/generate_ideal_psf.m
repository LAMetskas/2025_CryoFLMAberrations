function generate_ideal_psf(obj)
    % phaseretrieve_mle - Generate pupil function based on phase retrieval
    % algorithm described in the paper.
    % Edited by Lihongjia, 2024/05/07
    
    % Initialize variables
    z = obj.Zpos;
    k_z = compute_axial_frequency(obj);
    R = obj.PSFsize;
    
    

    % Finalize coefficients and calculate images
    [imgsOut_raw, pupil_phase, pupil_mag] = generate_psf(obj,obj.PSF_size_forMLE, z, k_z, R);
    
    save_results(obj, imgsOut_raw);
    %draw_and_save_results(obj, imgsOut_raw, obj.Mpsf_subroi);
    %draw_colormap(imgsOut_raw,obj.Mpsf_subroi);
    %    coeffVec_final'
    %draw_colormap_even_frames_together(imgsOut_raw,obj.Mpsf_subroi);

    %visualize_zx_plane(imgsOut_raw,obj.Mpsf_subroi);

    %visualize_yz_plane(imgsOut_raw,obj.Mpsf_subroi);
    % Update object properties
    %update_obj_properties(obj, coeffVec_final, pupil_phase, pupil_mag, imgsOut_raw);
end

function k_z = compute_axial_frequency(obj)
    % Compute axial frequency
    k_z = sqrt((obj.PRstruct.RefractiveIndex / obj.PRstruct.Lambda)^2 - obj.k_r.^2) .* obj.NA_constrain;
end


function [imgsOut_raw, pupil_phase, pupil_mag] = generate_psf(obj,R_measuredpsf, z, k_z, R)
    % Generate PSF images based on final coefficients
    [imgsOut_raw, aberrPhase, pupil_mag, pupil_phase] = calc_img_fromCoeffVec( z, obj.ZernikeorderN, k_z, obj.sVec, obj.bgVec, R, obj.Z.ZM, R_measuredpsf, obj.NA_constrain);
end


% function draw_and_save_results(obj, crop_fittedpsf_raw, crop_measuredpsf)
%     % Draw comparison and save results
%     draw_comparison(crop_fittedpsf_raw, 'fitted', crop_measuredpsf, 'measured');
%     draw_retrieved_psf(crop_fittedpsf_raw);
%     save2tiff(crop_fittedpsf_raw, 'final_fitted_psf_raw.tif');
%     save2tiff(crop_measuredpsf, 'final_measured_psf.tif');
%     save('no_aberration_psf.mat', 'crop_fittedpsf_raw');
%     save('measured_psf.mat', 'crop_measuredpsf');
% end
function save_results(obj, crop_fittedpsf_raw)
    % Draw comparison and save results
    %draw_comparison(crop_fittedpsf_raw, 'fitted', crop_measuredpsf, 'measured');
    %draw_retrieved_psf(crop_fittedpsf_raw);
    save2tiff(crop_fittedpsf_raw, 'final_fitted_psf_raw.tif');
    %save2tiff(crop_measuredpsf, 'final_measured_psf.tif');
    ims=crop_fittedpsf_raw;
    save('no_aberration_psf.mat', 'ims');
    %save('measured_psf.mat', 'crop_measuredpsf');
end

%function [imgsOut, aberrPhase, pupil_mag, pupil_phase] = calc_img_fromCoeffVec(coeffVec, zVec, ZernikeorderN, k_z, sVec, bgVec, gBlur, resizeFactor, R, ZM, PSF_size_forMLE, NA_constrain)
function [imgsOut, aberrPhase, pupil_mag, pupil_phase] = calc_img_fromCoeffVec(zVec, ZernikeorderN, k_z, sVec, bgVec, R, ZM, PSF_size_forMLE, NA_constrain)
    % Initialize variables
    N = length(zVec);
    Rpupil_mag = zeros(R, R, N); % retrieved pupil function magnitude
    Rpupil_phase = zeros(R, R, N); % retrieved pupil function phase
    pupil_mag = zeros(R, R);

    nZernikeCoeffs = (ZernikeorderN + 1) * (ZernikeorderN + 2) / 2; 
    
    % Expand the Zernike coefficients if necessary
    coeffVec_m = zeros(nZernikeCoeffs,1);
    coeffVec_m(1)=1;
    % Generate the pupil magnitude using Zernike polynomials
    for k = 1 : nZernikeCoeffs 
        pupil_mag = pupil_mag + ZM(:, :, k) .* coeffVec_m(k);
    end

    coeffVec_p=zeros(nZernikeCoeffs,1);
    coeffVec_p(1)=1;

    % Generate the pupil phase
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
            Iimg = Iimg * sVec(j) ;
            Iimg = Iimg + max(0, bgVec(j));
        end

        % Apply Gaussian blur if needed
        %if gBlur > 0
        %    h = fspecial('gaussian', [5 5], gBlur);
        %    Iimg = imfilter(Iimg, h);
        %end
        
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



function draw_retrieved_psf(imgsOut)
    % Create a full-screen figure
    h = figure('Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    
    % Number of frames in the PSF
    nFrames = size(imgsOut, 3);
    
    % Loop through and display each PSF image
    for i = 1:nFrames
        subplot(3, ceil(nFrames / 3), i);  % Position each image in a grid
        imshow(imgsOut(:,:,i), []);  % Display each image
        title(sprintf('Pos %d', i));  % Add a title to each subplot
    end
    
    % Save the figure as a .jpg file
    saveas(h, 'retrieved_psf.jpg', 'jpg');
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

function draw_colormap(stack1, stack2)
    % Normalize the stacks for better visualization (optional)
    stack1 = mat2gray(stack1);  % Normalize to [0,1]
    stack2 = mat2gray(stack2);  % Normalize to [0,1]
    
    
    % Set percentile values for intensity clipping (adjust these based on your data)
    lower_percentile = 1;  % Lower bound (1st percentile)
    upper_percentile = 99; % Upper bound (99th percentile)
    
    % Preallocate RGB stacks
    rgbStack1 = zeros(size(stack1,1), size(stack1,2), 3, size(stack1,3));
    rgbStack2 = zeros(size(stack2,1), size(stack2,2), 3, size(stack2,3));
    
    % Process stack1: Normalize each slice and clip intensity
    for i = 1:size(stack1, 3)
        slice = stack1(:,:,i);
        % Clip intensity values to a specific percentile range
        low_val = prctile(slice(:), lower_percentile);
        high_val = prctile(slice(:), upper_percentile);
        slice = min(max(slice, low_val), high_val);  % Clip to range
        slice = mat2gray(slice);  % Normalize to [0,1]
        rgbStack1(:,:,:,i) = ind2rgb(gray2ind(slice, 256), jet(256));  % Apply colormap
    end
    
    % Process stack2: Normalize each slice and clip intensity
    for i = 1:size(stack2, 3)
        slice = stack2(:,:,i);
        % Clip intensity values to a specific percentile range
        low_val = prctile(slice(:), lower_percentile);
        high_val = prctile(slice(:), upper_percentile);
        slice = min(max(slice, low_val), high_val);  % Clip to range
        slice = mat2gray(slice);  % Normalize to [0,1]
        rgbStack2(:,:,:,i) = ind2rgb(gray2ind(slice, 256), jet(256));  % Apply colormap
    end
    
    % Visualize all slices from stack1
    figure;
    montage(rgbStack1, 'Size', [1 NaN]);  % Display montage with 5 rows
    colorbar;  % Add colorbar for reference
    title('PSF Stack 1: All Slices with Clipping and Colormap');
    
    % Visualize all slices from stack2
    figure;
    montage(rgbStack2, 'Size', [1 NaN]);  % Display montage with 5 rows
    colorbar;  % Add colorbar for reference
    title('PSF Stack 2: All Slices with Clipping and Colormap');
end


function draw_colormap_even_frames_together(stack1, stack2)
    % Normalize the stacks for better visualization (optional)
    stack1 = mat2gray(stack1);  % Normalize to [0,1]
    stack2 = mat2gray(stack2);  % Normalize to [0,1]
    
    % Set percentile values for intensity clipping
    lower_percentile = 1;  % Lower bound (1st percentile)
    upper_percentile = 99; % Upper bound (99th percentile)
    
    % Preallocate RGB stacks for odd frames
    odd_indices = 2:2:size(stack1, 3);  % Indices of odd frames
    rgbStack1 = zeros(size(stack1,1), size(stack1,2), 3, length(odd_indices));  % For stack1
    rgbStack2 = zeros(size(stack2,1), size(stack2,2), 3, length(odd_indices));  % For stack2
    
    % Process stack1: Normalize each odd slice and clip intensity
    for idx = 1:length(odd_indices)
        i = odd_indices(idx);
        slice1 = stack1(:,:,i);
        
        % Clip intensity values to a specific percentile range for stack1
        low_val1 = prctile(slice1(:), lower_percentile);
        high_val1 = prctile(slice1(:), upper_percentile);
        slice1 = min(max(slice1, low_val1), high_val1);  % Clip to range
        slice1 = mat2gray(slice1);  % Normalize to [0,1]
        %rgbStack1(:,:,:,idx) = ind2rgb(gray2ind(slice1, 256), jet(256));  % Apply colormap
        rgbStack1(:,:,:,idx) = ind2rgb(gray2ind(slice1, 256), parula(256));  % Apply parula colormap
    end
    
    % Process stack2: Normalize each odd slice and clip intensity
    for idx = 1:length(odd_indices)
        i = odd_indices(idx);
        slice2 = stack2(:,:,i);
        
        % Clip intensity values to a specific percentile range for stack2
        low_val2 = prctile(slice2(:), lower_percentile);
        high_val2 = prctile(slice2(:), upper_percentile);
        slice2 = min(max(slice2, low_val2), high_val2);  % Clip to range
        slice2 = mat2gray(slice2);  % Normalize to [0,1]
        %rgbStack2(:,:,:,idx) = ind2rgb(gray2ind(slice2, 256), jet(256));  % Apply colormap
        rgbStack2(:,:,:,idx) = ind2rgb(gray2ind(slice2, 256), parula(256));  % Apply parula colormap
    end
    
    % Visualize all odd slices from stack1 in a montage
    figure;
    montage(rgbStack1, 'Size', [1 NaN]);  % Display the montage
    title('PSF Stack 1: Odd Frames');
    
    % Visualize all odd slices from stack2 in a montage
    figure;
    montage(rgbStack2, 'Size', [1 NaN]);  % Display the montage
    title('PSF Stack 2: Odd Frames');
end

function visualize_xz_plane(stack1, stack2)
    % Normalize the stacks for better visualization (optional)
    stack1 = mat2gray(stack1);  % Normalize to [0,1]
    stack2 = mat2gray(stack2);  % Normalize to [0,1]
    
    % Set percentile values for intensity clipping (adjust based on your data)
    lower_percentile = 1;  % Lower bound (1st percentile)
    upper_percentile = 99; % Upper bound (99th percentile)
    
    % Select a specific y-slice for the x-z plane visualization
    y_slice = round(size(stack1, 2) / 2);  % Take the middle y slice for x-z visualization

    % Extract the x-z planes for both stacks
    xz_plane_stack1 = squeeze(stack1(:, y_slice, :));  % Extract and squeeze to get 2D x-z image
    xz_plane_stack2 = squeeze(stack2(:, y_slice, :));  % Extract and squeeze to get 2D x-z image
    
    % Clip intensity values for xz_plane_stack1
    low_val1 = prctile(xz_plane_stack1(:), lower_percentile);
    high_val1 = prctile(xz_plane_stack1(:), upper_percentile);
    xz_plane_stack1 = min(max(xz_plane_stack1, low_val1), high_val1);  % Clip to range
    xz_plane_stack1 = mat2gray(xz_plane_stack1);  % Normalize to [0,1]
    rgb_xz_plane_stack1 = ind2rgb(gray2ind(xz_plane_stack1, 256), parula(256));  % Apply colormap

    % Clip intensity values for xz_plane_stack2
    low_val2 = prctile(xz_plane_stack2(:), lower_percentile);
    high_val2 = prctile(xz_plane_stack2(:), upper_percentile);
    xz_plane_stack2 = min(max(xz_plane_stack2, low_val2), high_val2);  % Clip to range
    xz_plane_stack2 = mat2gray(xz_plane_stack2);  % Normalize to [0,1]
    rgb_xz_plane_stack2 = ind2rgb(gray2ind(xz_plane_stack2, 256), parula(256));  % Apply colormap
    
    % Visualize the x-z plane for stack1
    figure;
    imshow(rgb_xz_plane_stack1);
    title('X-Z Plane of Stack 1');
    colorbar;

    % Visualize the x-z plane for stack2
    figure;
    imshow(rgb_xz_plane_stack2);
    title('X-Z Plane of Stack 2');
    colorbar;
end

function visualize_yz_plane(stack1, stack2)
    % Normalize the stacks for better visualization (optional)
    stack1 = mat2gray(stack1);  % Normalize to [0,1]
    stack2 = mat2gray(stack2);  % Normalize to [0,1]
    
    % Set percentile values for intensity clipping (adjust based on your data)
    lower_percentile = 1;  % Lower bound (1st percentile)
    upper_percentile = 99; % Upper bound (99th percentile)
    
    % Select a specific x-slice for the y-z plane visualization
    x_slice = round(size(stack1, 1) / 2);  % Take the middle x slice for y-z visualization

    % Extract the y-z planes for both stacks
    yz_plane_stack1 = squeeze(stack1(x_slice, :, :));  % Extract and squeeze to get 2D y-z image
    yz_plane_stack2 = squeeze(stack2(x_slice, :, :));  % Extract and squeeze to get 2D y-z image
    
    % Clip intensity values for yz_plane_stack1
    low_val1 = prctile(yz_plane_stack1(:), lower_percentile);
    high_val1 = prctile(yz_plane_stack1(:), upper_percentile);
    yz_plane_stack1 = min(max(yz_plane_stack1, low_val1), high_val1);  % Clip to range
    yz_plane_stack1 = mat2gray(yz_plane_stack1);  % Normalize to [0,1]
    rgb_yz_plane_stack1 = ind2rgb(gray2ind(yz_plane_stack1, 256), parula(256));  % Apply colormap

    % Clip intensity values for yz_plane_stack2
    low_val2 = prctile(yz_plane_stack2(:), lower_percentile);
    high_val2 = prctile(yz_plane_stack2(:), upper_percentile);
    yz_plane_stack2 = min(max(yz_plane_stack2, low_val2), high_val2);  % Clip to range
    yz_plane_stack2 = mat2gray(yz_plane_stack2);  % Normalize to [0,1]
    rgb_yz_plane_stack2 = ind2rgb(gray2ind(yz_plane_stack2, 256), parula(256));  % Apply colormap
    
    % Visualize the y-z plane for stack1
    figure;
    imshow(rgb_yz_plane_stack1);
    title('Y-Z Plane of Stack 1');
    colorbar;

    % Visualize the y-z plane for stack2
    figure;
    imshow(rgb_yz_plane_stack2);
    title('Y-Z Plane of Stack 2');
    colorbar;
end

function visualize_zx_plane(stack1, stack2)
%function visualize_xz_plane_with_z_as_vertical_larger(stack1, stack2)
    % Normalize the stacks for better visualization (optional)
    stack1 = mat2gray(stack1);  % Normalize to [0,1]
    stack2 = mat2gray(stack2);  % Normalize to [0,1]
    
    % Set percentile values for intensity clipping (adjust based on your data)
    lower_percentile = 1;  % Lower bound (1st percentile)
    upper_percentile = 99; % Upper bound (99th percentile)
    
    % Select a specific y-slice for the x-z plane visualization
    y_slice = round(size(stack1, 2) / 2);  % Take the middle y slice for x-z visualization

    % Extract the x-z planes for both stacks
    xz_plane_stack1 = squeeze(stack1(:, y_slice, :));  % Extract and squeeze to get 2D x-z image
    xz_plane_stack2 = squeeze(stack2(:, y_slice, :));  % Extract and squeeze to get 2D x-z image
    
    % Clip intensity values for xz_plane_stack1
    low_val1 = prctile(xz_plane_stack1(:), lower_percentile);
    high_val1 = prctile(xz_plane_stack1(:), upper_percentile);
    xz_plane_stack1 = min(max(xz_plane_stack1, low_val1), high_val1);  % Clip to range
    xz_plane_stack1 = mat2gray(xz_plane_stack1);  % Normalize to [0,1]
    rgb_xz_plane_stack1 = ind2rgb(gray2ind(xz_plane_stack1, 256), parula(256));  % Apply colormap

    % Clip intensity values for xz_plane_stack2
    low_val2 = prctile(xz_plane_stack2(:), lower_percentile);
    high_val2 = prctile(xz_plane_stack2(:), upper_percentile);
    xz_plane_stack2 = min(max(xz_plane_stack2, low_val2), high_val2);  % Clip to range
    xz_plane_stack2 = mat2gray(xz_plane_stack2);  % Normalize to [0,1]
    rgb_xz_plane_stack2 = ind2rgb(gray2ind(xz_plane_stack2, 256), parula(256));  % Apply colormap
    
    % Transpose the x-z plane to make Z vertical and X horizontal
    rgb_xz_plane_stack1 = permute(rgb_xz_plane_stack1, [2, 1, 3]);  % Reorder dimensions
    rgb_xz_plane_stack2 = permute(rgb_xz_plane_stack2, [2, 1, 3]);  % Reorder dimensions
    
    % Visualize the x-z plane for stack1 with Z as the vertical axis
    figure('Name', 'X-Z Plane Stack 1', 'Position', [100, 100, 1000, 800]);  % Create a large figure window
    imshow(rgb_xz_plane_stack1, 'InitialMagnification', 'fit');  % Adjust magnification
    title('X-Z Plane of Stack 1 (Z as vertical)');
    colorbar;

    % Visualize the x-z plane for stack2 with Z as the vertical axis
    figure('Name', 'X-Z Plane Stack 2', 'Position', [1200, 100, 1000, 800]);  % Create a large figure window
    imshow(rgb_xz_plane_stack2, 'InitialMagnification', 'fit');  % Adjust magnification
    title('X-Z Plane of Stack 2 (Z as vertical)');
    colorbar;
end



