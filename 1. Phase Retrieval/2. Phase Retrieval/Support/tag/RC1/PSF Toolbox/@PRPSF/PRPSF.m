classdef PRPSF < handle
    % PRPSF - Phase Retrieval PSF Model for Optical System Aberration Analysis
    %
    % Usage:
    %   obj = PRPSF();  % Create a new PRPSF object
    %
    % Description:
    %   This class supports MLE-based phase retrieval to extract pupil functions,
    %   Zernike phase coefficients, and generate PSFs for localization microscopy.
    %
    % Dependencies:
    %   - Zernike_Polynomials class
    %   - Image stacks (signal and noise)
    %   - FFT-based image processing

    %% ----------------------------
    % Public Properties (User Input)
    % -----------------------------
    properties
        % Z-stack range and indexing
        Zstart = -1.4;
        Zend = 1.4;
        Zstep = 0.2;
        Zindstart = 1;
        Zindend = 6;
        Zindstep = 1;
        DatadimX;
        DatadimY;
        DatadimZ

        % Data dimensions and raw data
        Beadcenter;         % Manually or automatically defined bead center
        BeadXYshift;        % Estimated XY shift from center
        BeadData;           % 3D signal image stack (photon converted)
        NoiseData;          % 3D background image stack (photon converted)

        % Imaging parameters
        Pixelsize;          % Microns per pixel
        PSFsize;            % Size of output PSF
        SubroiSize;         % Size of region around bead
        OTFratioSize;       % Size used in OTF estimation
        ZernikeorderN;      % Order of Zernike decomposition
        PSF_size_forMLE;    % Size of PSF used in MLE cost function
        resizeFactor;       % FFT upsampling factor
        gBlur;              % Gaussian blur sigma
        CCDoffset;          % Camera baseline offset (ADU)
        Gain;               % Electron-to-ADU gain

        % File and directory control
        FileDir;            % Path to original data
        FileName;           % Filename
        SaveDir;            % Where to store output
        
        % Phase retrieval setup
        IterationNum;       % Total iteration count
        IterationNumK;      % Iteration count after which data preprocessing changes
        log_index = 0;      % Index for logging internal metrics

        % MLE/Phase Retrieval Inputs
        PRstruct;           % Structure for PR results
        PSFstruct;          % Structure for generated PSFs
        Zpos;               % Z positions in microns
        coeffVecStart;      % Initial Zernike coefficients
        bgVec;              % Estimated background (photons/frame)
        sVec;               % Estimated signal (photons/frame)
        r;                  % Pupil radius (in pixels)
        Z;                  % Zernike_Polynomials object

        % Zernike coefficients (phase and magnitude)
        phase_coefficients;
        mag_coefficients;
        phase_coefficients_pupil;
    end

    %% ----------------------------
    % Private Internal Properties
    % ----------------------------
    properties (SetAccess = private, GetAccess = private)
        PhiC;       % Angular map for raw data
        ZoC;        % Radial map for raw data
        Phi;        % Angular map for PSF
        k_r;        % Normalized spatial frequency
        Zo;         % Radial map for PSF
        NA_constrain; % Fourier support mask based on NA
    end

    %% ----------------------------
    % Accessible Intermediate Results
    % ----------------------------
    properties (SetAccess = public, GetAccess = public)
        Mpsf_subroi;
        Mpsf_subroi_input;
        Mpsf_extend;
        ncc_values;
        mse_values;
        loglikelihood_values;
    end

    %% ----------------------------
    % Constructor
    % ----------------------------
    %% ----------------------------
    % Methods
    % ----------------------------
    methods
        function obj = PRPSF()
            % PRPSF Constructor
        end

        function prepdata(obj)
            % Convert ADU to photons and extract valid bead signal
            in = (obj.BeadData - obj.CCDoffset) ./ obj.Gain;
            [obj.DatadimY, obj.DatadimX, obj.DatadimZ] = size(in);
            obj.Mpsf_subroi_input = in;
            obj.Mpsf_extend = in;
        end

        function precomputeParam(obj)
            % Generate spatial maps for Fourier optics and Zernike decomposition
            [XC, YC] = meshgrid(-obj.DatadimX/2:obj.DatadimX/2-1, -obj.DatadimY/2:obj.DatadimY/2-1);
            obj.PhiC = atan2(YC, XC);
            obj.ZoC = hypot(XC, YC);

            [X, Y] = meshgrid(-obj.PSFsize/2:obj.PSFsize/2-1);
            obj.Zo = hypot(X, Y);
            scale = obj.PSFsize * obj.Pixelsize;
            obj.k_r = obj.Zo / scale;
            obj.Phi = atan2(Y, X);

            Freq_max = obj.PRstruct.NA / obj.PRstruct.Lambda;
            obj.NA_constrain = obj.k_r < Freq_max;

            zk = Zernike_Polynomials();
            zk.Ordering = 'Noll';
            zk.setN(obj.ZernikeorderN);
            zk.initialize();
            [Zrho, Ztheta, Zinit] = zk.params3_Zernike(obj.Phi, obj.k_r, obj.PRstruct.NA, obj.PRstruct.Lambda);
            zk.matrix_Z(Zrho, Ztheta, Zinit);
            obj.Z = zk;
        end

        function datapreprocess(obj)
            % Crop and center PSF using estimated XY shift
            realsize0 = floor(obj.SubroiSize / 2);
            realsize1 = ceil(obj.SubroiSize / 2);
            MpsfC = zeros(obj.SubroiSize, obj.SubroiSize, obj.DatadimZ);
            shiftxy = obj.Beadcenter - obj.BeadXYshift;

            for ii = 1:obj.DatadimZ
                tmp = fftshift(ifft2(obj.Mpsf_subroi_input(:,:,ii)));
                shiftphase = -obj.ZoC / obj.DatadimX .* cos(obj.PhiC) .* (obj.DatadimX/2 - shiftxy(1) - 1) ...
                             - obj.ZoC / obj.DatadimY .* sin(obj.PhiC) .* (obj.DatadimY/2 - shiftxy(2) - 1);
                tmp1 = fft2(tmp .* exp(-2 * pi * shiftphase * 1i));
                MpsfC(:,:,ii) = abs(tmp1(obj.DatadimY/2-realsize0+1:obj.DatadimY/2+realsize1, ...
                                          obj.DatadimX/2-realsize0+1:obj.DatadimX/2+realsize1));
            end
            obj.Mpsf_subroi = MpsfC;
        end

        function genMpsf(obj)
            % Crop centered subregion to MLE-compatible size
            mid = floor(size(obj.Mpsf_subroi,1)/2);
            half_size = obj.PSF_size_forMLE / 2;
            idx = mid-half_size+1 : mid+half_size;
            obj.Mpsf_subroi = obj.Mpsf_subroi(idx, idx, :);
        end

        function estdata(obj)
            % Estimate signal and background vectors from input stacks
            disp('Calculating background...');
            bgData = (obj.NoiseData - obj.CCDoffset) ./ obj.Gain;
            obj.bgVec = mean(bgData, [1, 2]);

            disp('Calculating signal...');
            r_psf = size(obj.Mpsf_subroi, 1);
            signalData = obj.Mpsf_subroi - obj.bgVec .* ones(r_psf, r_psf, obj.DatadimZ);
            obj.sVec = sum(reshape(signalData, [], obj.DatadimZ));

            % disp(mean(obj.bgVec));
            % disp(mean(obj.sVec));
        end

        function saveObj(obj, saveName)
            % Save the PRPSF object to a .mat file
            if isempty(obj.SaveDir)
                error('Save directory not specified.');
            end
            PRFileName = fullfile(obj.SaveDir, [obj.FileName(1:end-4), saveName, '.mat']);
            save(PRFileName, 'obj');
        end

        function plot_similarity(obj, data, title_n, t)
            % Plot similarity metrics (NCC, MSE, log-likelihood, etc.)
            figure;
            plot(data, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
            title(t);
            xlabel('Iteration');
            ylabel('Metric Value');
            grid on;
            saveas(gcf, [title_n, '.png']);
            close(gcf);
        end
    end
end