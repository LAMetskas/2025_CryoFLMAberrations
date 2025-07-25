% Generate a PSF using Zernike approximation.
function PSF = Generate_PSF(PRT)

   Ordering = 'Wyant';
   %Ordering = 'Noll';
   %fprintf('ordering = %s\n', Ordering);

   CN_complex = PRT.CN_complex;
   lambda = PRT.Lambda;
   Maglateral = PRT.Magnify;
   n = PRT.RefractiveIndex;
   NA = PRT.NA;
   pixel_size = PRT.PixelSizeCCD;
   R = PRT.OutputPsfSize;
   Zernike_nmax = PRT.ZernikeOrder;
   Zstack = PRT.Zstack;

   % Precompute the Zernike polynomial coefficients.
   Z = Zernike_Polynomials();
   Z.Ordering = Ordering;
%  if strcmp(Ordering, 'Wyant')
%     Z.Wnmax = Zernike_nmax;
%     Z.initialize_Wyant();
%  else   % Noll
%     %Zernike_lmax = Z.n_Zernikes_Noll(Zernike_nmax);
%     %Z.Nlmax = Zernike_lmax;
%     Z.Nnmax = Zernike_nmax;
%     Z.initialize_Noll();
%  end
   Z.setN(Zernike_nmax);
   Z.initialize();

   % Compute parameters associated with the Zernike polynomials.
   [rho, theta, NA_constraint, k_r] = ...
      Z.params4_Zernike(R, pixel_size, Maglateral, NA, lambda);

   % Generate the complex pupil function using the given parameters for
   % (rho, theta, init, w) where w are the weights and init the nonzero pixels.
%  if strcmp(Ordering, 'Wyant')
%     complex_MagZ = Z.poly_sum_Wyant(rho, theta, NA_constraint, CN_complex);
%  else
%     complex_MagZ = Z.poly_sum_Noll(rho, theta, NA_constraint, CN_complex);
%  end
   complex_MagZ = Z.poly_sum(rho, theta, NA_constraint, CN_complex);

   % Generate PSF from the generated pupil function.
   N = length(Zstack);
   PSF = zeros(R, R, N);
   k_z = sqrt((n/lambda)^2 - k_r.^2) .* NA_constraint;
   C = 2 .* pi .* 1i .* k_z;
   for j = 1 : N
      %defocus_phase = 2 * pi * Zstack(j) .* k_z .* 1i;
      defocus_phase = C * Zstack(j);
      pupil_complex = complex_MagZ .* exp(defocus_phase);
      psfA = abs(fftshift(fft2(pupil_complex)));
      Fig2 = psfA.^2;
      PSF(:, :, j) = Fig2 ./ R^4;
   end

end
