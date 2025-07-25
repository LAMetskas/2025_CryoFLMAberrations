% This is an example script for using the PRPSF class

%% select a measured PSF data
% test data : PSF Toolbox\test data\darkredbead-2015-10-27-11-6-28.mat
[FileName, FileDir] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
F = load(fullfile(FileDir,FileName));
sz=size(F.dataset);

%% create object and set input properties of PRPSF class
% here the settings are for SEQ microscope
probj=PRPSF();
probj.CCDoffset=100;
probj.Gain=2;
probj.PRstruct.NA=1.49;
probj.PRstruct.Lambda=0.69;
probj.PRstruct.RefractiveIndex=1.51;
probj.Pixelsize=0.1;
probj.PSFsize=128;
probj.SubroiSize=50;
probj.OTFratioSize=60;
probj.ZernikeorderN=7;
probj.Zstart = F.Params.Zrange(1); %micron
probj.Zend = F.Params.Zrange(3); %micron
probj.Zstep = F.Params.Zrange(2); %micron
probj.Zindstart = 1; %index
probj.Zindend = sz(4);
probj.Zindstep = 1;
probj.IterationNum=25;
probj.IterationNumK=5;

%% generate PR result
probj.FileDir=FileDir;
probj.FileName=FileName;
probj.prepdata();
probj.precomputeParam();
probj.findXYshift();
probj.datapreprocess();
probj.genMpsf();
probj.phaseretrieve();
probj.genZKresult();
probj.findOTFparam();
% generate figures for phase retrieval results
probj.genPRfigs('PSF');
probj.genPRfigs('pupil');
probj.genPRfigs('zernike');
probj.calcrlb();
%% save the PRPSF object
resdir=fullfile(FileDir,'PR_result\');
if ~exist(resdir,'dir')
    mkdir(resdir)
end
probj.SaveDir=resdir;
probj.saveObj('_PR_result_test');