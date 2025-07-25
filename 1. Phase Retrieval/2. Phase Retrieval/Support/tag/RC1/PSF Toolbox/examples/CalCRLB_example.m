% This is an example script for using the CalCRLB class

addpath('C:\Users\xufbi\Documents\GitHub\tag\RC1\PSF Toolbox');

%% method 1 ------------------calculate CRLB from phase retrieval result ----------
% test data: PSF Toolbox\test data\darkredbead-2015-10-27-11-6-28_PR_result_test.mat
[FileName, FileDir] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
% load the selected PRPSF object
tmp=load(fullfile(FileDir,FileName));
namei=fields(tmp);
probj=tmp.(namei{1});

% create PSF_pupil object and set input properties of PSF_pupil class
crobj=CalCRLB(probj.PRstruct,'pupil');
Num=20;
crobj.Pixelsize=0.1;%micron
crobj.Xpos=0 + 0.*rand(Num,1);
crobj.Ypos=0 + 0.*rand(Num,1);
crobj.Zpos=linspace(-1,1,Num)';
crobj.Photon=1000.*ones(Num,1);
crobj.Bg=20.*ones(Num,1);
crobj.Boxsize=16;
crobj.Deltax=0.1;%pixel
crobj.Deltaz=0.01;%micron

crobj.prepInputparam();
crobj.calcrlb();
crobj.genfigs();

%% method 2 ------------------calculate CRLB with user defined PSF model ----------
% generate a PRstruct if only use the zernike coefficients for PSF generation
% in this case, create object by: crobj=CalCRLB(PRstruct,'zernike')
PRstruct.NA=1.49;
PRstruct.Lambda=0.69;
PRstruct.RefractiveIndex=1.52;
PRstruct.Pupil.phase=zeros(128,128);
PRstruct.Pupil.mag=zeros(128,128);
PRstruct.Zernike_phase=[0,0,0,0,1,0,0,0,0];% this generate astigmatism aberration
PRstruct.Zernike_mag=[1,0,0,0,0,0,0,0,0];
PRstruct.SigmaX=2;
PRstruct.SigmaY=2;

crobj=CalCRLB(PRstruct,'zernike');
Num=20;
crobj.Pixelsize=0.1;%micron
crobj.Xpos=0 + 0.*rand(Num,1);
crobj.Ypos=0 + 0.*rand(Num,1);
crobj.Zpos=linspace(-0.5,0,Num)';
crobj.Photon=1000.*ones(Num,1);
crobj.Bg=20.*ones(Num,1);
crobj.Boxsize=16;
crobj.Deltax=0.1;%pixel
crobj.Deltaz=0.01;%micron
crobj.PSFobj.PSFsize=128;
crobj.PSFobj.nMed = 1.33;

crobj.prepInputparam();
crobj.calcrlb();
crobj.genfigs();
