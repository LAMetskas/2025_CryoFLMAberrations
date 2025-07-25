% This is an example script for using the OptimPR_Bi class

%% select a measured PSF data
% test data : PSF Toolbox\test data\darkredbead-2015-8-11-11-23-33.mat
[FileName, FileDir] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
F = load(fullfile(FileDir,FileName));
sz=size(F.dataset);
%% create object and set input properties of PRPSF class
% here the settings are for TIRF microscope
oprobj=OptimPR_Bi();
oprobj.PRobjA.PRobj.CCDoffset=106.2;
oprobj.PRobjA.PRobj.Gain=10.41;
oprobj.PRobjA.PRobj.PRstruct.NA=1.46;
oprobj.PRobjA.PRobj.PRstruct.Lambda=0.69;
oprobj.PRobjA.PRobj.PRstruct.RefractiveIndex=1.52;
oprobj.PRobjA.PRobj.Pixelsize=0.1067; % micron
oprobj.PRobjA.PRobj.PSFsize=128;
oprobj.PRobjA.PRobj.SubroiSize=60;
oprobj.PRobjA.PRobj.OTFratioSize=60;
oprobj.PRobjA.PRobj.ZernikeorderN=7;
oprobj.PRobjA.PRobj.Zstart = F.Params.Zrange(1); %micron
oprobj.PRobjA.PRobj.Zend = F.Params.Zrange(3); %micron
oprobj.PRobjA.PRobj.Zstep = F.Params.Zrange(2); %micron
oprobj.PRobjA.PRobj.Zindstart = 1; %index
oprobj.PRobjA.PRobj.Zindend = sz(4);
oprobj.PRobjA.PRobj.Zindstep = 1;
oprobj.PRobjA.PRobj.IterationNum=25;
oprobj.PRobjA.PRobj.IterationNumK=5;
oprobj.PRobjA.PRobj.PRstruct.SigmaX=2;
oprobj.PRobjA.PRobj.PRstruct.SigmaY=2;
oprobj.PRobjA.FileDir=FileDir;
oprobj.PRobjA.FileName=FileName;
oprobj.PlaneDis=0.43;% um
oprobj.IterationMonte=50;
oprobj.FileDir=FileDir;
oprobj.FileName=FileName;

%% generate initial PR result
oprobj.prepdata('EMCCD');% 'sCMOS' or 'EMCCD'
oprobj.initialPR();
% show PR result from plane A
oprobj.PRobjA.PRobj.genPRfigs('zernike');
oprobj.PRobjA.PRobj.genPRfigs('pupil');
oprobj.PRobjA.PRobj.genPRfigs('PSF');
oprobj.PRobjA.PRobj.calcrlb();
% show PR result from plane B
oprobj.PRobjB.PRobj.genPRfigs('zernike');
oprobj.PRobjB.PRobj.genPRfigs('pupil');
oprobj.PRobjB.PRobj.genPRfigs('PSF');
oprobj.PRobjB.PRobj.calcrlb();

%% optimization of PR result
oprobj.PRobjA.FitZrange=[oprobj.PRobjA.PRobj.Zstart,oprobj.PRobjA.PRobj.Zend-oprobj.PlaneDis];
oprobj.optimPR();
oprobj.PRobjA.genPR();
oprobj.PRobjA.genOTFpsf();
oprobj.PRobjB.genPR();
oprobj.PRobjB.genOTFpsf();

% show optimized PR result
oprobj.MCResult.LLTrace
oprobj.PRobjA.PRobj.PRstruct
oprobj.PRobjB.PRobj.PRstruct
    % PR result from plane A
oprobj.PRobjA.PRobj.genPRfigs('zernike');
oprobj.PRobjA.PRobj.genPRfigs('pupil');
oprobj.PRobjA.PRobj.genPRfigs('PSF');
oprobj.PRobjA.PRobj.calcrlb();
    % PR result from plane B
oprobj.PRobjB.PRobj.genPRfigs('zernike');
oprobj.PRobjB.PRobj.genPRfigs('pupil');
oprobj.PRobjB.PRobj.genPRfigs('PSF');
oprobj.PRobjB.PRobj.calcrlb();

%% generate parameters for 3D localization
oprobj.BoxSizeFit=16;
oprobj.CenterSize=4;
oprobj.findZestimator();

%% PR fit bead data
oprobj.PRobjA.genSamplePSF();
oprobj.FitType=1;%1: for 'EMCCD' and 2: for 'sCMOS', however, 'sCMOS' fitting is not available yet
oprobj.Iterationsfit=100;
oprobj.fitPSFdata();
% show fitting result
oprobj.fitFigure;

%% save OptimPR_Bi object
resdir=fullfile(FileDir,'PR_result\');
if ~exist(resdir,'dir')
    mkdir(resdir)
end
oprobj.SaveDir=resdir;
oprobj.SaveName='_PR_result_optimBi_test';
oprobj.saveObj();

