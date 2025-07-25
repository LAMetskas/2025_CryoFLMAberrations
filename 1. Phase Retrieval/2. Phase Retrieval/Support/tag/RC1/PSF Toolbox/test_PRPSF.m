me=userpath;
cd([me(1:end-1),'\SM_Analysis\development\PSF Toolbox\'])


%% test PRPSF

probj=PRPSF();
probj.CCDoffset=130;
probj.Gain=0.075;
probj.PRstruct.NA=1.49;
probj.PRstruct.Lambda=0.69;
probj.PRstruct.RefractiveIndex=1.51;
probj.Pixelsize=0.1;
probj.PSFsize=128;
probj.SubroiSize=50;
probj.OTFratioSize=60;
probj.ZernikeorderN=7;
probj.Zstart = -1; %micron
probj.Zend = 1; %micron
probj.Zstep = 0.1; %micron
probj.Zindstart = 1; %index
probj.Zindend = 21;
probj.Zindstep = 2;
probj.IterationNum=25;
probj.IterationNumK=5;

[FileName FileDir filterInd] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
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

probj.genPRfigs('PSF');
probj.genPRfigs('pupil');
probj.genPRfigs('zernike');
probj.calcrlb();
%%
resdir=fullfile(FileDir,'PR_result\');
if ~exist(resdir,'dir')
    mkdir(resdir)
end
probj.SaveDir=resdir;
probj.saveObj('_PR_result');
%% test PSF_pupil
psfobj=PSF_pupil(probj.PRstruct);
Num=20;
psfobj.Pixelsize=0.1;
psfobj.precomputeParam();

psfobj.Xpos=-2 + 0.*rand(Num,1);
psfobj.Ypos=-2 + 0.*rand(Num,1);
psfobj.Zpos=-0.5 + 0.*rand(Num,1);
psfobj.Boxsize=32;
psfobj.genPSF();
psfobj.scalePSF();

dipshow(psfobj.PSFs)
dipshow(psfobj.ScaledPSFs)

%% test CalCRLB
crobj=CalCRLB(probj.PRstruct,'pupil');
Num=20;
crobj.Pixelsize=0.1;%micron
crobj.Xpos=0 + 0.*rand(Num,1);
crobj.Ypos=0 + 0.*rand(Num,1);
crobj.Zpos=linspace(-1,1,Num)';
crobj.Photon=1000.*ones(Num,1);
crobj.Bg=40.*ones(Num,1);
crobj.Boxsize=16;
crobj.Deltax=0.1;%pixel
crobj.Deltaz=0.01;%micron

crobj.prepInputparam();
crobj.calcrlb();
crobj.genfigs();

%% test optimPR_Ast
oprobj=OptimPR_Ast();

[FileName FileDir filterInd] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
F = load(fullfile(FileDir,FileName));
GainFileDir='Y:\sCMOS Calibrations\RB mic\Slow readout mode';
GainFileName='GainCalibration-2015-9-3-16-15-4.mat';

oprobj.PRobj.CCDoffset=130;
oprobj.PRobj.Gain=0.075;
oprobj.PRobj.PRstruct.NA=1.49;
oprobj.PRobj.PRstruct.Lambda=0.69;
oprobj.PRobj.PRstruct.RefractiveIndex=1.51;
oprobj.PRobj.Pixelsize=0.096;
oprobj.PRobj.PSFsize=128;
oprobj.PRobj.SubroiSize=40;
oprobj.PRobj.OTFratioSize=60;
oprobj.PRobj.ZernikeorderN=7;
oprobj.PRobj.Zstart = F.Params.Zrange(1); %micron
oprobj.PRobj.Zend = F.Params.Zrange(end); %micron
oprobj.PRobj.Zstep = 0.1; %micron
oprobj.PRobj.Zindstart = 1; %index
oprobj.PRobj.Zindend = numel(F.Params.Zrange);
oprobj.PRobj.Zindstep = 4;
oprobj.PRobj.IterationNum=25;
oprobj.PRobj.IterationNumK=5;

oprobj.FileDir=FileDir;
oprobj.FileName=FileName;
oprobj.GainFileDir=GainFileDir;
oprobj.GainFileName=GainFileName;

oprobj.prepdata('sCMOS');
oprobj.initialPR();
oprobj.FitZrange=[0.2,0.6];
oprobj.PRobj.PRstruct.SigmaX=2;
oprobj.PRobj.PRstruct.SigmaY=2;
oprobj.IterationMonte=50;
oprobj.optimPR();
oprobj.genPR();
oprobj.genOTFpsf();
oprobj.MCResult.LLTrace
oprobj.PRobj.PRstruct

oprobj.PRobj.genPRfigs('zernike');
oprobj.PRobj.genPRfigs('pupil');
oprobj.PRobj.genPRfigs('PSF');
oprobj.PRobj.calcrlb();
% PR fit bead data
oprobj.BoxSizeFit=16;% boxsize must be 16
oprobj.genSamplePSF();
oprobj.Iterationsfit=100;
[P,LL,SSE,crlb,CG]=oprobj.fitPSFdata();

oprobj.fitFigure(P,LL);









