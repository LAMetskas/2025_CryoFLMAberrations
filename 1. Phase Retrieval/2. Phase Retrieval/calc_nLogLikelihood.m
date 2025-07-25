% Copyright (c)2017, The Board of Trustees of The Leland Stanford Junior University.
% All rights reserved.
% Licensed under BSD 3-Clause License.

function [O] = calc_nLogLikelihood(obj, data, model)
% CALC_NLOGLIKELIHOOD - Compute the negative log-likelihood (NLL) under a Poisson noise model.
% Used in phase retrieval optimization.
%
% Inputs:
%   obj   - Object with logging and PSF metadata
%   data  - Measured PSF image stack (HxWxD)
%   model - Simulated PSF image stack (HxWxD)
%
% Output:
%   O     - Scalar value of negative log-likelihood

    % Compute NLL for Poisson model: sum(model - data * log(model))
    O = sum(model(:) - data(:) .* log(model(:)));

    % Periodically log NCC and MSE (every 100 iterations)
    if obj.log_index ~= 0 && mod(obj.log_index, 100) == 0
        ncc_values = calculate_ncc_for_stacks(data, model);
        ncc = mean(ncc_values);
        mse_err = immse(data, model);
        obj.ncc_values = [obj.ncc_values, ncc];
        obj.mse_values = [obj.mse_values, mse_err];
        obj.loglikelihood_values = [obj.loglikelihood_values, O];
    end

    % Increment log index
    obj.log_index = obj.log_index + 1;
end

function ncc_values = calculate_ncc_for_stacks(stack1, stack2)
% CALCULATE_NCC_FOR_STACKS - Compute NCC values slice-by-slice between two 3D stacks
% Inputs:
%   stack1, stack2 - Image stacks of same size (HxWxD)
% Output:
%   ncc_values     - 1xD vector of NCCs per slice

    if ~isequal(size(stack1), size(stack2))
        error('Input stacks must have the same size.');
    end

    num_images = size(stack1, 3);
    ncc_values = zeros(1, num_images);

    for i = 1:num_images
        ncc_values(i) = calculate_ncc(stack1(:, :, i), stack2(:, :, i));
    end
end

function ncc_value = calculate_ncc(img1, img2)
% CALCULATE_NCC - Compute normalized cross-correlation between two 2D images
% Inputs:
%   img1, img2 - 2D image slices
% Output:
%   ncc_value  - Scalar NCC

    img1 = double(img1);
    img2 = double(img2);

    img1_centered = img1 - mean(img1(:));
    img2_centered = img2 - mean(img2(:));

    numerator = sum(img1_centered(:) .* img2_centered(:));
    denominator = sqrt(sum(img1_centered(:).^2) * sum(img2_centered(:).^2));

    ncc_value = numerator / denominator;
end