% This is an example script for using the PSF_pupil class

%% select a saved phase retrieval result file
% test data: PSF Toolbox\test data\darkredbead-2015-10-27-11-6-28_PR_result_test.mat
[FileName, FileDir] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
% load the selected PRPSF object
tmp=load(fullfile(FileDir,FileName));
namei=fields(tmp);
probj=tmp.(namei{1});
%% create PSF_pupil object and set input properties of PSF_pupil class
psfobj=PSF_pupil(probj.PRstruct);
% generate 20 PSF images with given x, y, and z positions
Num=20;
psfobj.Xpos=-1 + 2.*rand(Num,1);% pixel
psfobj.Ypos=-1 + 2.*rand(Num,1);% pixel
psfobj.Zpos=-0.5 + 0.5.*rand(Num,1);% micron
psfobj.Boxsize=32;
psfobj.Pixelsize=0.1; % micron
psfobj.nMed=1.33;

psfobj.precomputeParam();
psfobj.genPSF();
psfobj.scalePSF();
% show generated PSF images
dipshow(psfobj.PSFs)
dipshow(psfobj.ScaledPSFs)

%% generate PSFs considering index mismatch aberration
Num=20;
psfobj.Xpos=0.*rand(Num,1);% pixel
psfobj.Ypos=0.*rand(Num,1);% pixel
psfobj.Zpos=3;% micron
psfobj.ZposMed=linspace(-1,1,Num);% micron
psfobj.Boxsize=32;
psfobj.Pixelsize=0.1; % micron
psfobj.nMed=1.33;

psfobj.precomputeParam();
psfobj.genIMMPSF();
% show generated PSF images
dipshow(psfobj.IMMPSFs)