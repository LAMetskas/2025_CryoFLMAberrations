classdef PRPSF < handle
    % PRPSF class for generating phase retrieved pupil function and PSF
    %   Create object: obj = PRPSF();
    %
    % PRPSF Properties (Input):
    %   CCDoffset - Camera offset value(s)
    %   FileDir - Directory of the measured PSF data
    %   FileName - Name of the measured PSF data file
    %   Gain - Gain value(s) from camera calibration
    %   IterationNum - Number of iterations for phase retrieval algorithm
    %   IterationNumK - Number of iterations after which data preprocessing will change
    %   OTFratioSize - Size of the subregion for OTF rescale
    %   PSFsize - Output size of PSF
    %   Pixelsize - Pixel size at sample plane (microns)
    %   SubroiSize - Size of subregion for phase retrieval
    %   SaveDir - Directory to save the PRPSF object
    %   ZernikeorderN - Maximum order of Zernike coefficient
    %   Zstart, Zend, Zstep - Z range and step size (microns)
    %   Zindstart, Zindend, Zindstep - Indices for phase retrieval Z positions
    %
    % PRPSF Properties (Output):
    %   PRstruct - Structure containing phase retrieval results
    %   PSFstruct - Structure containing various PSFs
    %
    % PRPSF Methods:
    %   prepdata - Convert ADU count to photon count and average over time
    %   precomputeParam - Generate images for k-space operation
    %   datapreprocess - Shift the bead to the center and crop image to subregions
    %   genMpsf - Generate normalized measured PSF for phase retrieval
    %   saveObj - Save PRPSF object
    %   calcrlb - Calculate CRLB based on PSF model
    %   phaseretrieve_mle - Generate pupil function using phase retrieval
    %   findOTFparam - Find SigmaX and SigmaY for OTF rescale
    %   findXYshift - Find x,y shift of selected bead image
    %   fitdefocus - Find amount of defocus in z of the measured PSF
    %   genPRfigs - Generate figures of phase retrieval result
    %   genZKresult - Expand phase retrieved pupil function into Zernike polynomials

    properties
        % Input properties
        Zstart = -1; % Start position of z, in microns
        Zend = 1; % End position of z, in microns
        Zstep = 0.4; % Step size of z, in microns
        Zindstart = 1; % Start index for phase retrieval
        Zindend = 6; % End index for phase retrieval
        Zindstep = 1; % Index step for phase retrieval
        Beadcenter; % Pixel position of selected bead 
        BeadXYshift; % Difference between Beadcenter and found bead position, in pixels
        BeadData; % Bead images after ADU-to-photon conversion and time averaging
        DatadimX; % X dimension of BeadData
        DatadimY; % Y dimension of BeadData
        DatadimZ; % Z dimension of BeadData
        PSFsize; % Output size of PSF
        SubroiSize; % Subregion size for phase retrieval
        Pixelsize; % Pixel size at sample plane, in microns
        IterationNum; % Iteration number for phase retrieval algorithm
        IterationNumK; % Iteration number for PR algorithm after data preprocessing change 
        ZernikeorderN; % Maximum order of Zernike coefficient (Wyant ordering)
        OTFratioSize; % Subregion size for OTF rescale
        PRstruct; % Structure for phase retrieval results
        PSFstruct; % Structure for various PSFs
        CCDoffset; % Camera offset value(s)
        Gain; % Gain value(s) from camera calibration
        FileDir; % Directory of measured PSF data
        FileName; % Name of the measured PSF data file
        SaveDir; % Save directory for PRPSF object
        Z; % Zernike_Polynomials object
        Enableunwrap; % 1: Unwrap phase before Zernike expansion, 0: Do not unwrap phase
        Zpos; % Z positions for images (added by Fan Xu)
        coeffVecStart; % Initial coefficients for Zernike polynomials (MLE algorithm)
        bgVec; % Background vector for the bead image stack (in photons)
        sVec; % Signal vector for the bead image stack (in photons)
        NoiseData; % Noise images after ADU-to-photon conversion and time averaging
        r; % Radius of pupil in pixels (MLE algorithm)
        gBlur; % Sigma for Gaussian filter (MLE algorithm)
        resizeFactor; % Upsample scale for FFT (MLE algorithm)
        PSF_size_forMLE; % PSF size for MLE algorithm
        mag_coefficients; % Magnitude coefficients for Zernike expansion
        phase_coefficients; % Phase coefficients for Zernike expansion
    end

    properties (SetAccess = private, GetAccess = private)
        % Precomputed images for k-space operation
        PhiC; % Phi coordinates for BeadData, size DatadimX x DatadimY
        ZoC; % R coordinates for BeadData, size DatadimX x DatadimY
        Phi; % Phi coordinates for output PSF, size PSFsize x PSFsize
        k_r; % k_r coordinates for output OTF, size PSFsize x PSFsize
        Zo; % R coordinates for output PSF, size PSFsize x PSFsize
        NA_constrain; % Circular function defining cut-off frequency in Fourier space
    end
    
    % Fan Xu change 'SetAccess' to public 
    properties (SetAccess = public, GetAccess = public)
        Mpsf_subroi; % Measured PSF after data preprocessing, 3D image stack SubroiSize x SubroiSize x DatadimZ
        Mpsf_extend; % Normalized measured PSF for phase retrieval, 3D image stack PSFsize x PSFsize x DatadimZ       
        log_index;
        ncc_values;
        mse_values;
        loglikelihood_values;
    end
    
    methods
        function obj = PRPSF()
            % Constructor for PRPSF class
        end

        function prepdata(obj)
            % prepdata - Convert ADU count to photon count and average over time dimension
            in = (obj.BeadData - obj.CCDoffset) .* obj.Gain;
            [obj.DatadimY, obj.DatadimX, obj.DatadimZ] = size(in);
            obj.Mpsf_subroi = in;

            %prepared data for XY shift estimation
            R=size(obj.Mpsf_subroi,1);
            N=size(obj.Mpsf_subroi,3);
            MpsfL=zeros(R,R,N);
            for ii=1:N
                Mpsfo=obj.Mpsf_subroi(:,:,ii);  %change
                %Edge=[mean(Mpsfo(1:3,:)),mean(Mpsfo(R-3:R,:)),mean(Mpsfo(:,1:3)),mean(Mpsfo(:,R-3:R))];
                Edge=[mean(Mpsfo(1,:)),mean(Mpsfo(R,:)),mean(Mpsfo(:,1)),mean(Mpsfo(:,R))];
                bg=max(Edge);
                Fig2=(Mpsfo-bg);
                Fig2(Fig2<=0)=0;
                MpsfL(:,:,ii)=Fig2;
            end

            obj.Mpsf_extend=MpsfL;
            %obj.Mpsf_extend=obj.Mpsf_subroi;
        end

        function precomputeParam(obj)
            % precomputeParam - Generate images for k-space operation and save in precomputed parameters.
            [XC, YC] = meshgrid(-obj.DatadimX/2:obj.DatadimX/2-1, -obj.DatadimY/2:obj.DatadimY/2-1);
            obj.PhiC = atan2(YC, XC);
            obj.ZoC = hypot(XC, YC); % Use hypot for better numerical stability

            [X, Y] = meshgrid(-obj.PSFsize/2:obj.PSFsize/2-1, -obj.PSFsize/2:obj.PSFsize/2-1);
            obj.Zo = hypot(X, Y);
            scale = obj.PSFsize * obj.Pixelsize;
            obj.k_r = obj.Zo / scale;
            obj.Phi = atan2(Y, X);
            
            Freq_max = obj.PRstruct.NA / obj.PRstruct.Lambda;
            obj.NA_constrain = obj.k_r < Freq_max;

            % Create Zernike_Polynomials object
            zk = Zernike_Polynomials();
            zk.Ordering = 'Noll';
            %zk.Ordering = 'Wyant';
            zk.setN(obj.ZernikeorderN);
            zk.initialize();
            [Zrho, Ztheta, Zinit] = zk.params3_Zernike(obj.Phi, obj.k_r, obj.PRstruct.NA, obj.PRstruct.Lambda);
            zk.matrix_Z(Zrho, Ztheta, Zinit);
            obj.Z = zk;
        end

        function datapreprocess(obj)
            % datapreprocess - Shift bead to the center and crop image to subregions.
            realsize0 = floor(obj.SubroiSize / 2);
            realsize1 = ceil(obj.SubroiSize / 2);
            MpsfC = zeros(obj.SubroiSize, obj.SubroiSize, obj.DatadimZ);
            shiftxy = obj.Beadcenter - obj.BeadXYshift;

            for ii = 1:obj.DatadimZ
                tmp = fftshift(ifft2(squeeze(obj.Mpsf_subroi(:, :, ii))));
                shiftphase = -obj.ZoC / obj.DatadimX .* cos(obj.PhiC) .* (obj.DatadimX / 2 - shiftxy(1) - 1) ...
                             - obj.ZoC / obj.DatadimY .* sin(obj.PhiC) .* (obj.DatadimY / 2 - shiftxy(2) - 1);
                tmp1 = fft2(tmp .* exp(-2 * pi * shiftphase * 1i));
                startx = -realsize0 + obj.DatadimX / 2 + 1;
                endx = realsize1 + obj.DatadimX / 2;
                starty = -realsize0 + obj.DatadimY / 2 + 1;
                endy = realsize1 + obj.DatadimY / 2;
                MpsfC(:, :, ii) = abs(tmp1(starty:endy, startx:endx));
            end

            obj.Mpsf_subroi = MpsfC;
        end
        
        %for loglikelihood
        function genMpsf(obj)
            % genMpsf - Generate normalized measured PSF for phase retrieval
            %obj.Mpsf_extend = obj.Mpsf_subroi;
            
            %crop data based on the MLE_size
            Measured_PSF_size = size(obj.Mpsf_subroi,1);
            obj.Mpsf_subroi = obj.Mpsf_subroi(Measured_PSF_size/2-obj.PSF_size_forMLE/2+1 : Measured_PSF_size/2+obj.PSF_size_forMLE/2, ...
                                 Measured_PSF_size/2-obj.PSF_size_forMLE/2+1 : Measured_PSF_size/2+obj.PSF_size_forMLE/2, :);
        end
        
        %for mse
        function genMpsf_bak2(obj)
            % genMpsf- operate on Mpsf_subroi, generate normalized measured PSF that are used for phase
            % retrieval. 
            R=size(obj.Mpsf_subroi,1);
            N=size(obj.Mpsf_subroi,3);
            MpsfL=zeros(R,R,N);
            for ii=1:N
                Mpsfo=obj.Mpsf_subroi(:,:,ii);  %change
                %Edge=[mean(Mpsfo(1:3,:)),mean(Mpsfo(R-3:R,:)),mean(Mpsfo(:,1:3)),mean(Mpsfo(:,R-3:R))];
                Edge=[mean(Mpsfo(1,:)),mean(Mpsfo(R,:)),mean(Mpsfo(:,1)),mean(Mpsfo(:,R))];
                %Fig2=(Mpsfo-obj.bgVec(:,ii));
                %Fig2=Mpsfo;
                bg=max(Edge);
                Fig2=(Mpsfo-bg);
                Fig2(Fig2<=0)=0;
                MpsfL(:,:,ii)=Fig2;
            end

            obj.Mpsf_extend=MpsfL;
        end

        function estdata(obj)
            % estdata - Calculate the total photons and background for each frame
            nFrames = size(obj.NoiseData, 3);
            disp('Calculating background...');           

            % Convert from counts to photons
            bgData = (obj.NoiseData - obj.CCDoffset) .* obj.Gain;

            %preproce
            %bgData_=obj.remove_poisson_noise_stack(bgData);
            obj.bgVec = mean(bgData, [1, 2]);

            disp('... Background calculated.');

            disp('Calculating signal...');
            r_psf = size(obj.Mpsf_subroi, 1);

            %obj.Mpsf_subroi=obj.remove_poisson_noise_stack(obj.Mpsf_subroi);
            %signalData_ = obj.Mpsf_subroi - obj.bgVec .* ones(r_psf, r_psf, nFrames);
            %obj.sVec = sum(reshape(signalData_, [], nFrames));
            signalData = obj.Mpsf_subroi - obj.bgVec .* ones(r_psf, r_psf, nFrames);
            obj.sVec = sum(reshape(signalData, [], nFrames));
            
            %disp(sum(obj.bgVec))
            %disp(sum(obj.sVec))
            disp('... Signal calculated.');
        end

        function saveObj(obj, saveName)
            % saveObj - Save PRPSF object
            filename = obj.FileName;
            if isempty(obj.SaveDir)
                error('PRPSF:NoSaveDir','Save directory is empty')
            else
                PRFileName = fullfile(obj.SaveDir, [filename(1:end-4) saveName '.mat']);
                save(PRFileName, 'obj');
            end
        end

        function generate_pm(obj)
            % generate_pm - Generate phase and magnitude pupil functions
            R = 256;
            %n = 11; % Replace with appropriate logic if n should be determined dynamically
            n=size(obj.phase_coefficients,2);
            n=ceil((obj.ZernikeorderN+1)*(obj.ZernikeorderN+2)/2)
            pupil_phasenorm = zeros(R, R);

            for k = 1:n
                pupil_phasenorm = pupil_phasenorm + obj.Z.ZM(:, :, k) .* obj.phase_coefficients(k);
            end

            pupil_magnorm = zeros(R, R);
            obj.mag_coefficients = zeros(n);
            obj.mag_coefficients(1) = 0;

            for k = 1:n
                pupil_magnorm = pupil_magnorm + obj.Z.ZM(:, :, k) .* obj.mag_coefficients(k);
            end

            obj.PRstruct.Fittedpupil.mag = pupil_magnorm;
            obj.PRstruct.Fittedpupil.phase = pupil_phasenorm;
        end

        function calcrlb(obj)
            % calcrlb - Calculate CRLB based on PSF model from phase retrieval result
            crobj = CalCRLB(obj.PRstruct, 'pupil');
            z = obj.Zpos;
            Num = numel(z);
            crobj.Pixelsize = obj.Pixelsize; % microns
            crobj.Xpos = zeros(Num, 1);
            crobj.Ypos = zeros(Num, 1);
            crobj.Zpos = z';
            crobj.Photon = 1000 .* ones(Num, 1);
            crobj.Bg = 2 .* ones(Num, 1);
            crobj.Boxsize = 16;
            crobj.Deltax = 0.1; % pixel
            crobj.Deltaz = 0.01; % micron
            
            crobj.prepInputparam();
            crobj.calcrlb();
            crobj.genfigs();
        end

        function plot_similarity(obj, data, title_n)
            % plot_similarity - Plot similarity metrics (e.g., NCC, MSE)
            figure;
            plot(data, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
            title(title_n);
            xlabel('Index');
            ylabel('Value');
            grid on;
            axis tight;
            %saveas(gcf, fullfile('./Intermediate results/', [title_n, '.png']));
            saveas(gcf, [title_n, '.png']);
            close(gcf);
        end
        function transformed_stack = remove_poisson_noise_stack(obj,image_stack)
            % remove_poisson_noise_stack - Remove Poisson noise from an image stack
            %
            % Inputs:
            %   image_stack - The input image stack with Poisson noise (3D array)
            %
            % Outputs:
            %   denoised_stack - The denoised image stack with reduced Poisson noise (3D array)
        
            % Get the number of slices in the stack
            num_slices = size(image_stack, 3);
            
            % Initialize the denoised stack
            transformed_stack = zeros(size(image_stack));
        
            % Loop through each slice in the stack
            for i = 1:num_slices
                % Extract the i-th slice
                image_slice = image_stack(:, :, i);
        
                % Apply Anscombe transformation
                transformed_slice = obj.anscombe_transform(image_slice);
        
                transformed_stack(:, :, i) = transformed_slice;
            end
        
            % Display the first slice of the results for quick visualization
            %figure;
            %subplot(1, 2, 1); imshow(image_stack(:, :, 1), []); title('Noisy Image (First Slice)');
            %subplot(1, 2, 2); imshow(denoised_stack(:, :, 1), []); title('Denoised Image (First Slice)');
        end
        function transformed_image = anscombe_transform(obj, image)
            % anscombe_transform - Apply the Anscombe transformation to stabilize variance
            %
            % Inputs:
            %   image - The input image with Poisson noise
            %
            % Outputs:
            %   transformed_image - The image after Anscombe transformation
        
            transformed_image = 2 * sqrt(image + 3/8);
        end

    end
end
