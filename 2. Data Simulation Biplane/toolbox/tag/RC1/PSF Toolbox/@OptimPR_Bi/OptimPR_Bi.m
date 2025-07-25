classdef OptimPR_Bi < handle
    % OptimPR_Bi class for optimizing PR result which will be used for 3D
    % localization based on dual focal plane method
    %   create object: obj = OptimPR_Bi();
    %
    % OptimPR_Bi Properties (Implemented object):
    %   OTFobj - 
    %   PRobjA - 
    %   PRobjB - 
    %   PSFobj - 
    %
    % OptimPR_Bi Properties (Input):
    %   CenterSize - 
    %   FileDir - 
    %   FileName - 
    %   FitType - 
    %   FitZrange - 
    %   GainFileDir - 
    %   GainFileName - 
    %   IterationMonte - 
    %   Iterationsfit - 
    %   SaveDir - 
    %   SaveName - 
    %   
    % OptimPR_Bi Properties (Output):
    %   Fitresult - 
    %   Iratio - 
    %   MCResult - 
    %   PlaneDis - 
    %   PlaneShift - 
    %   Zestimator - 
    %
    % OptimPR_Bi Methods:
    %   prepdata - converting ADU count to photon count and averaging over time dimension
    %   initialPR - generate initial phase retrieval result
    %   optimPR - optimize phase retrieval result
    %   findZestimator - find coefficients for initial estimation of z positions in 3D localization
    %   fitFigure - generate figures from localization results of measured PSFs
    %   fitPSFdata - find 3D localization results of the measured PSFs that were used for phase retrieval
    %   objectivefun - the objective function to be optimized in the optimization step
    %   runMCMC - generate optimized parameters for phase retrieval using Monte Carlo Markov Chain       
    %   saveObj - save OptimPR_Ast object in SaveDir with a SaveName
    %
    %   see also OTFrescale OptimPR_Ast PSF_pupil
    properties
        BoxSizeFit;% subregion size for 3D localization, it must be 16
        CenterSize;% subregion size used for estimation of x,y,z positions of the measured PSFs
        FileDir;% file directory of the measured PSF data
        FileName;% file name of the measured PSF data
        FitType;% fitting type for specific camera, options are 1: for 'EMCCD' and 2: for 'sCMOS'
        GainFileDir;% file directory of gain calibration result
        GainFileName;% file name of gain calibration result
        IterationMonte;% number of iterations in optimization of phase retrieved PSF
        Iterationsfit;% number of iterations in 3D localization
        Iratio;% photon count ratio between the PSFs measured in plane A and plane B
        % MCResult - a structure of optimization history
        %   LLTrace: log likelihood values from accepted step
        MCResult;
        % Fitresult - fitting result of bead data that were used for phase retrieval
        %   P: results of fitting parameters, which is [x; y; I; Bg; z]
        %   LL: log likelihood between data and PSF model
        %   SSE: sum of square error between data and PSF model
        %   CRLB: CRLB of each fitting parameters
        %   CG: convergence of each fitting parameters
        Fitresult;
        OTFobj;% object of OTFrescale class
        PlaneDis;% plane separation in z between plane A and plane B, unit is micron
        PlaneShift;% plane shift in xy between plane A and plane B, unit is pixel
        PRobjA;% object of OptimPR_Ast class, used for phase retrieval of measured PSFs from plane A 
        PRobjB;% object of OptimPR_Ast class, used for phase retrieval of measured PSFs from plane B
        PSFobj;% object of PSF_pupil class
        SaveDir;% save directory of OptimPR_Bi object
        SaveName;% save name of OptimPR_Bi object
        Zestimator;% coefficients from polynomial fitting for initial estimation of z positions in 3D localization
    end
    
    methods
        function obj=OptimPR_Bi()
            obj.PRobjA=OptimPR_Ast;
            obj.PRobjB=OptimPR_Ast;
            obj.OTFobj=OTFrescale();
        end
        function prepdata(obj,cameratype)
            % prepdata - converting ADU count to photon count and averaging
            % over time dimension. 
            %   It is depend on the camera type because of the difference in
            %   gain calibration results. It calls the 'prepdata' method of
            %   PRobjA, then make a copy of all properties in PRobjA to
            %   PRobjB.
            %   Input parameter: cameratype - type of camera, options are
            %   'EMCCD' and 'sCMOS'
            %   
            %   see also OptimPR_Ast.prepdata
            obj.PRobjA.prepdata(cameratype);
            fieldn=fieldnames(obj.PRobjA.PRobj);
            for ii=1:numel(fieldn)
                try
                obj.PRobjB.PRobj.(fieldn{ii})=obj.PRobjA.PRobj.(fieldn{ii});
                catch me
                end
            end
        end
        function initialPR(obj)
            % initialPR - generate initial phase retrieval result and use
            % the resulting zernike coefficient to minimize the x,y and z
            % shifts. 
            %   It uses the method function 'initialPR' in
            %   OptimPR_Ast class for both plane A and plane B, in which,
            %   plane B has a initial z shift relative to plane A. Then, the
            %   x,y,z shift between two planes are calculated from the
            %   initial phase retrieval result
            %
            %   see also OptimPR_Ast.initialPR

           obj.PRobjA.initialPR();
           obj.PSFobj=PSF_pupil(obj.PRobjA.PRobj.PRstruct);
           obj.PRobjB.PRobj.Zstart=obj.PRobjA.PRobj.Zstart+obj.PlaneDis;
           obj.PRobjB.PRobj.Zend=obj.PRobjA.PRobj.Zend+obj.PlaneDis;
           obj.PRobjB.initialPR();
           obj.PlaneDis=obj.PRobjB.PRobj.Zstart-obj.PRobjA.PRobj.Zstart;
           z=[obj.PRobjA.PRobj.Zstart:obj.PRobjA.PRobj.Zstep:obj.PRobjA.PRobj.Zend];
           [val,indf]=min(abs(z+obj.PlaneDis/2));
           img1=obj.PRobjA.PRobj.Mpsf_subroi;
           img2=obj.PRobjB.PRobj.Mpsf_subroi;
           c=findshift(img2(:,:,indf),img1(:,:,indf),'iter'); %position(img2)-position(img1)
           obj.PlaneShift=c; 
           IratioL=squeeze(sum(sum(img2,1),2))./squeeze(sum(sum(img1,1),2));
           obj.Iratio=mean(IratioL); 
        end
        
        function optimPR(obj)
            % optimPR - optimize phase retrieval result. 
            %   It uses Monte Carlo Markov Chain to iteratively update 5
            %   parameters, which are the refractive index, sigmax and sigmay
            %   in OTF rescale, 'SubroiSize' and 'IterationNumK' defined in
            %   PRPSF class, in order to minimize the difference between
            %   phase retrieved PSFs and measured PSFs.
            %
            %   see also runMCMC objectivefun
            [est,dLL]=obj.runMCMC();
            fieldn={'PRobjA','PRobjB'};
            for ii=1:numel(fieldn)
                obj.(fieldn{ii}).PRobj.PRstruct.RefractiveIndex=est(3);
                obj.(fieldn{ii}).PRobj.PRstruct.SigmaX=est(4);
                obj.(fieldn{ii}).PRobj.PRstruct.SigmaY=est(5);
                obj.(fieldn{ii}).PRobj.SubroiSize=est(1);
                obj.(fieldn{ii}).PRobj.IterationNumK=est(2);
            end
            if isfield(obj.MCResult,'LLTrace')
                obj.MCResult.LLTrace=cat(1,obj.MCResult.LLTrace,dLL);
            else
                obj.MCResult.LLTrace=dLL;
            end
        end
        
       function FileOut=saveObj(obj)
            filename=obj.FileName;
            
            if isempty(obj.SaveDir)
                error('PR:NoSaveDir','save directry is empty')
            else
                PRFileName=fullfile(obj.SaveDir,[filename(1:end-4) obj.SaveName '.mat']);
                save(PRFileName,'obj');
            end
            
            FileOut=[filename(1:end-4) obj.SaveName '.mat'];
        end
        
    end
    
end

