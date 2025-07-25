function findXYshift(obj)
    % findXYshift - Determine x, y shift of the selected bead image using a 2D Gaussian
    % fit to the most in-focus PSF. The selected bead image is cropped around the
    % selected pixel from the BeadData.
    
    R = [obj.DatadimX, obj.DatadimY]; % Image dimensions
    midZ = ceil(obj.DatadimZ / 2);    % Midpoint in the Z-dimension
    
    % Extract the central slice of the 3D PSF
    Fig1 = squeeze(obj.Mpsf_extend(:, :, midZ));
    
    % Display the image and get user-selected coordinates
    h = dipshow(Fig1);
    maxdim = max(size(Fig1));
    diptruesize(h, 1000 / maxdim * 100);  % Adjust image size for display
    Centers = dipgetcoords(1);  % Get the coordinates from user selection
    close(h)
    
    Ri = 32;  % Crop size for FFT-based operations
    PHI = obj.PhiC;  % Predefined angle parameter
    Zo = obj.ZoC;    % Predefined distance parameter
    
    % Define cropping indices for FFT operation
    realsize0 = floor(Ri / 2);
    realsize1 = ceil(Ri / 2);
    starty = -realsize0 + R(2) / 2 + 1;
    endy = realsize1 + R(2) / 2;
    startx = -realsize0 + R(1) / 2 + 1;
    endx = realsize1 + R(1) / 2;
    
    % Compute the initial FFT and apply phase shift correction
    tmp = fftshift(ifft2(Fig1));
    shiftphase = -Zo ./ R(1) .* cos(PHI) .* (R(1) / 2 - Centers(1,1) - 1) - ...
                 Zo ./ R(2) .* sin(PHI) .* (R(2) / 2 - Centers(1,2) - 1);
    tmp1 = fft2(tmp .* exp(-2 * pi .* shiftphase .* 1i));
    
    % Crop the FFT result to focus on the region of interest
    Mfocus = abs(tmp1(starty:endy, startx:endx));
    
    % Perform 2D Gaussian fitting using fminsearch
    [Dxy, ~, ~] = fminsearch(@(x) gaussianD(x, Mfocus, Ri), [1000, 1, 1, 0.5, 0.5]);
    
    % Store the results in the object properties
    obj.Beadcenter = Centers;
    obj.BeadXYshift = Dxy(4:5);
    %obj.BeadXYshift = [0, 0];  % This line seems to reset the shift; consider revising if needed
end

function [sse, data] = gaussianD(x, input_im, R)
    % gaussianD - Calculate the sum of squared errors (sse) between the input image
    % and a 2D Gaussian model. This function is used for optimization.
    
    I = x(1);        % Intensity scaling factor
    sigma = x(2);    % Gaussian standard deviation
    bg = x(3);       % Background level
    x0 = x(4);       % X-shift
    y0 = x(5);       % Y-shift
    
    % Generate a grid for the Gaussian model
    [xx, yy] = meshgrid(-R/2:R/2-1, -R/2:R/2-1);
    
    % Compute the 2D Gaussian model
    Model = I .* exp(-xx.^2 ./ (2 * sigma^2)) .* exp(-yy.^2 ./ (2 * sigma^2)) + bg;
    
    % Calculate phase shift correction
    PHI = atan2(yy, xx);
    Zo = sqrt(xx.^2 + yy.^2);
    tmp = fftshift(ifft2(input_im));
    shiftphase = -Zo ./ R .* cos(PHI) .* x0 - Zo ./ R .* sin(PHI) .* y0;
    tmp1 = fft2(tmp .* exp(-2 * pi .* shiftphase .* 1i));
    
    % Extract the magnitude of the FFT result
    data = abs(tmp1);
    
    % Calculate sum of squared errors between model and data
    sse = sum(sum((Model - data).^2));
end

