% This is an example script for using the PSF_zernike class

%% method 1 ------------------generate PSF from phase retrieval result ----------
% select a saved phase retrieval result file
% test data: PSF Toolbox\test data\darkredbead-2015-10-27-11-6-28_PR_result_test.mat
[FileName, FileDir] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
% load the selected PRPSF object
tmp=load(fullfile(FileDir,FileName));
namei=fields(tmp);
probj=tmp.(namei{1});
%% create PSF_pupil object and set input properties of PSF_pupil class
psfobj=PSF_zernike(probj.PRstruct);
% generate 20 PSF images with given x, y, and z positions
Num=20;
psfobj.Xpos=-2 + 2.*rand(Num,1);% pixel
psfobj.Ypos=-2 + 2.*rand(Num,1);% pixel
psfobj.Zpos=-0.5 + 0.5.*rand(Num,1);% micron
psfobj.Boxsize=32;
psfobj.PSFsize=128;
psfobj.Pixelsize=0.1; % micron
psfobj.nMed=1.33;

psfobj.precomputeParam();
psfobj.genPupil();
psfobj.genPSF();
psfobj.scalePSF();
% show generated PSF images
dipshow(psfobj.PSFs)
dipshow(psfobj.ScaledPSFs)
% show pupil image
dipshow(psfobj.Pupil.mag)
dipshow(psfobj.Pupil.phase)

%% method 2 ------------------generate PSF with user defined PSF model ----------
% generate a PRstruct if only use the zernike coefficients for PSF generation
% in this case, create object by: psfobj=PSF_zernike(PRstruct)
PRstruct.NA=1.49;
PRstruct.Lambda=0.69;
PRstruct.RefractiveIndex=1.52;
PRstruct.Pupil.phase=zeros(128,128);
PRstruct.Pupil.mag=zeros(128,128);
PRstruct.Zernike_phase=[0,0,0,0,1,0,0,0,0];% this generate astigmatism aberration
PRstruct.Zernike_mag=[1,0,0,0,0,0,0,0,0];
PRstruct.SigmaX=2;
PRstruct.SigmaY=2;

psfobj=PSF_zernike(PRstruct);
% generate 20 PSF images with given x, y, and z positions
Num=20;
psfobj.Xpos=-2 + 2.*rand(Num,1);% pixel
psfobj.Ypos=-2 + 2.*rand(Num,1);% pixel
psfobj.Zpos=-0.5 + 0.5.*rand(Num,1);% micron
psfobj.Boxsize=32;
psfobj.Pixelsize=0.1; % micron
psfobj.PSFsize=128;
psfobj.nMed=1.33;

psfobj.precomputeParam();
psfobj.genPupil();
tic
psfobj.genPSF();
psfobj.scalePSF();
toc
% show generated PSF images
dipshow(psfobj.PSFs)
dipshow(psfobj.ScaledPSFs)
% show pupil image
dipshow(psfobj.Pupil.mag)
dipshow(psfobj.Pupil.phase)

%% generate PSFs considering index mismatch aberration, with user defined PSF model 
PRstruct.NA=1.49;
PRstruct.Lambda=0.69;
PRstruct.RefractiveIndex=1.52;
PRstruct.Pupil.phase=zeros(128,128);
PRstruct.Pupil.mag=zeros(128,128);
PRstruct.Zernike_phase=[0,0,0,0,1,0,0,0,0];% this generate astigmatism aberration
PRstruct.Zernike_mag=[1,0,0,0,0,0,0,0,0];
PRstruct.SigmaX=2;
PRstruct.SigmaY=2;

psfobj=PSF_zernike(PRstruct);

Num=20;
psfobj.Xpos=0.*rand(Num,1);% pixel
psfobj.Ypos=0.*rand(Num,1);% pixel
psfobj.Zpos=3;% micron
psfobj.ZposMed=linspace(-1,1,Num);% micron
psfobj.Boxsize=32;
psfobj.Pixelsize=0.1; % micron
psfobj.PSFsize=128;
psfobj.nMed=1.33;

psfobj.precomputeParam();
psfobj.genPupil();
psfobj.genIMMPSF();
psfobj.scalePSF();
% show generated PSF images
dipshow(psfobj.IMMPSFs)


