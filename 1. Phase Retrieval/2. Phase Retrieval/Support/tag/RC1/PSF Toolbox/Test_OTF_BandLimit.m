
close all
clear all
clear classes
%%
a=OTF_BandLimit; % define the object
% set the properties
a.Bpp=0.096; % micron (the only input we want to rely on is the Bpp and the Data we get)
% PRobj Data:
load('Y:\Farzin\Sequential Imaging Data\2015_7_6_CRLBbeads\PR_result\darkredbead-2015-7-6-16-5-46_PR_optim_single.mat')

a.RawData=obj.Data(:,:,1); % left: defined object above.name of the property in class / right: obj from PR results.Data
a.RawData=mean(obj.Data,3); % left: defined object above.name of the property in class / right: obj from PR results.Data

a.findOTF_Radius() % finds radius of OTF of the data and then gives back NA_data
a

%% Synthetic PSF
M=[1,1];
Magnify=250; 
pixelsizeCCD=24;% micrion
boxsize=128;
n_med=1.33;
n_imm=1.52;
Intensity=2000;
BG1=0;
Refindex=[n_imm,n_med];
lambda=0.69;% micrion
NA=1.4946;
OTFparam1=[1,100,100,0]; % [magnitude,sigmaX,sigmaY,bg]
R=128;
% depth~=0, then if our model is correct we wouldn't have symmetry above and under zero in our synthetic PSF so it can work as a test of the model.
depth=single(0.0001); % micrion  
Zstack=[-1:0.1:1]';
[PSF,aber_phase]=IndexMismatchPSF(Refindex,NA,lambda,R,boxsize,pixelsizeCCD,...
                    Magnify,Zstack,depth,OTFparam1,Intensity,BG1);
            
PSF = single(noise(10000*PSF,'poisson'));

dipshow(PSF)
a.RawData=PSF(:,:,1);
%a.RawData=mean(PSF,3);

a.findOTF_Radius()
a

%%

OTF = ft(PSF(:,:,11))
figure
LRM = log(radialmean(abs(OTF)))

