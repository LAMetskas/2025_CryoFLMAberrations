function denoised_image = remove_poisson_noise(obj,image)
    % remove_poisson_noise - Remove Poisson noise from an image
    %
    % Inputs:
    %   image - The input image with Poisson noise
    %
    % Outputs:
    %   denoised_image - The denoised image with reduced Poisson noise

    % Apply Anscombe transformation
    transformed_image = anscombe_transform(image);
    
    % Apply denoising in the transformed domain
    denoised_transformed_image = denoise_image(transformed_image);
    
    % Apply inverse Anscombe transformation
    denoised_image = inverse_anscombe_transform(denoised_transformed_image);

    % Display the results
    figure;
    subplot(1, 2, 1); imshow(image, []); title('Noisy Image');
    subplot(1, 2, 2); imshow(denoised_image, []); title('Denoised Image');
end

function transformed_image = anscombe_transform(image)
    % anscombe_transform - Apply the Anscombe transformation to stabilize variance
    %
    % Inputs:
    %   image - The input image with Poisson noise
    %
    % Outputs:
    %   transformed_image - The image after Anscombe transformation

    transformed_image = 2 * sqrt(image + 3/8);
end

function denoised_image = denoise_image(image)
    % denoise_image - Denoise the image using a Gaussian filter
    %
    % Inputs:
    %   image - The input image to be denoised
    %
    % Outputs:
    %   denoised_image - The denoised image

    % Apply a Gaussian filter with a standard deviation (sigma) of 1
    sigma = 1;
    denoised_image = imgaussfilt(image, sigma);
end


function inv_transformed_image = inverse_anscombe_transform(image)
    % inverse_anscombe_transform - Apply the inverse Anscombe transformation
    %
    % Inputs:
    %   image - The input image after denoising (Anscombe domain)
    %
    % Outputs:
    %   inv_transformed_image - The image after inverse Anscombe transformation

    inv_transformed_image = (image / 2).^2 - 1/8;
end
