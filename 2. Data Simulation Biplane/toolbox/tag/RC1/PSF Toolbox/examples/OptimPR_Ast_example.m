% This is an example script for using the OptimPR_Ast class

%% select a measured PSF data
% test data : PSF Toolbox\test data\darkredbead-2015-10-27-11-6-28.mat
[FileName, FileDir] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
F = load(fullfile(FileDir,FileName));
sz=size(F.dataset);
% gainfile: PSF Toolbox\test data\GainCalibration-2015-8-17-18-7-15.mat
[GainFileName,GainFileDir] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
%% create object and set input properties of PRPSF class
% here the settings are for SEQ microscope
oprobj=OptimPR_Ast();
oprobj.PRobj.CCDoffset=100;
oprobj.PRobj.Gain=2;
oprobj.PRobj.PRstruct.NA=1.49;
oprobj.PRobj.PRstruct.Lambda=0.69;
oprobj.PRobj.PRstruct.RefractiveIndex=1.52;
oprobj.PRobj.Pixelsize=0.1;
oprobj.PRobj.PSFsize=128;
oprobj.PRobj.SubroiSize=50; % we cut the initial image. Different sizes might give better or worse results.
oprobj.PRobj.OTFratioSize=60;
oprobj.PRobj.ZernikeorderN=7;
oprobj.PRobj.Zstart = F.Params.Zrange(1); %micron
oprobj.PRobj.Zend = F.Params.Zrange(3); %micron
oprobj.PRobj.Zstep = F.Params.Zrange(2); %micron
oprobj.PRobj.Zindstart = 1; %index
oprobj.PRobj.Zindend = sz(4);
oprobj.PRobj.Zindstep = 1;
oprobj.PRobj.IterationNum=25;
oprobj.PRobj.IterationNumK=5;
oprobj.FitZrange=[-0.8,0.8];
oprobj.PRobj.PRstruct.SigmaX=2;
oprobj.PRobj.PRstruct.SigmaY=2;
oprobj.IterationMonte=50;

oprobj.FileDir=FileDir;
oprobj.FileName=FileName;
oprobj.GainFileDir=GainFileDir;
oprobj.GainFileName=GainFileName;

%% generate initial PR result
oprobj.prepdata('sCMOS');
oprobj.initialPR();
oprobj.PRobj.genPRfigs('zernike');
oprobj.PRobj.genPRfigs('pupil');
oprobj.PRobj.genPRfigs('PSF');
oprobj.PRobj.calcrlb();

%% optimization of PR result
oprobj.optimPR();
oprobj.genPR();
oprobj.genOTFpsf();
% show optimized PR result
oprobj.MCResult.LLTrace
oprobj.PRobj.PRstruct

oprobj.PRobj.genPRfigs('zernike');
oprobj.PRobj.genPRfigs('pupil');
oprobj.PRobj.genPRfigs('PSF');
oprobj.PRobj.calcrlb();

%% generate parameters for 3D localization
oprobj.FitZrange=[-0.8,0.8]; % change z range if necessary
oprobj.findZestimator();
oprobj.findAstparam();

%% fit bead data with PR PSF
oprobj.genSamplePSF();
oprobj.FitType=2;%1: for 'EMCCD' and 2: for 'sCMOS'
oprobj.Iterationsfit=100;
oprobj.fitPSFdata();
% show fitting result
oprobj.fitFigure();

%% save OptimPR_Ast object
resdir=fullfile(FileDir,'PR_result\');
if ~exist(resdir,'dir')
    mkdir(resdir)
end
oprobj.SaveDir=resdir;
oprobj.SaveName='_PR_result_optimAst_test';
oprobj.saveObj();






