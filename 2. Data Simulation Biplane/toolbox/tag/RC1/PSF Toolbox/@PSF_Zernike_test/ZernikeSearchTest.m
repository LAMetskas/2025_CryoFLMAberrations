% Direct search for matching data PSF with Zernike terms: 
%(fminsearch doesn't work good! Let's use MonteCarlo method)

%To use fminsearch we need initial value: take zernike coeffs of your PR data (based on Adaptive Optics, they should be enough)
%and use them as initial value:
clear all
clc
addpath('C:\Users\Farzin\Dropbox\Iso TIRF SR\Farzin Code')
NA=1.49;lambda=0.69;n_med=1.33;n_imm=1.52;OTFparam1=[1,2,2,0];R=128; 
Zstack=[-1:0.4:1];
%1- if we just use the initial values of our own / fourth term is focus
Theta0=[1.49,0.637,n_med,n_imm,OTFparam1(1),0.2,0,0,0,0,zeros(1,76)]; % we should get 41 inputs
zernike_psf=Zernike_PSF(Theta0,Zstack); 
dipshow(zernike_psf)

%2- if we use the initial values from the Phase Retrieval result on the Raw Data:
addpath('Y:\Farzin\Sequential Imaging Data\2015_7_30_CRLBpaper_test2_Collar19end\PR_result')
load('darkredbead2_1-2015-7-30-16-48-9_PR_optim_single.mat')
RawData=obj.Data; % 128-128-21
error=zeros(R,R,length(Zstack));
for ii=1:length(Zstack)
ScaledRawData(:,:,ii)=RawData(:,:,ii)/sum(sum(RawData(:,:,ii)));
end
ZernikePhaseCoeffs=obj.Results.PhaseCoeff;
ZernikeAmplitudeCoeffs=obj.Results.MagCoeff;
Opts = optimset()
% Theta0 from PhaseRetrieval:
Theta0=[[NA,lambda,n_med,n_imm,OTFparam1(1)],[ZernikePhaseCoeffs]];
zernike_psf=Zernike_PSF(Theta0,Zstack); 
dipshow(zernike_psf)
%Theta=[NA,lambda,n_med,n_imm,OTFparam1(1),Coeff_phase];
ErrFunc = @(Theta0,ScaledRawData,Zstack) mse(Zernike_PSF(Theta0,Zstack),ScaledRawData); % here the PSF is our zernike_psf
Opts = optimset() %Set options
% we use fminsearch to minimize the difference between RealData and Model
[Theta_Out,FVAL,EXITFLAG,OUTPUT] = fminsearch(ErrFunc,Theta0,Opts,ScaledRawData,Zstack) % Our initial guess:
Model = Zernike_PSF(Theta_Out,Zstack) % rewrite the initial zernike_psf with the optimized zernike phase coeffs
joinchannels('RGB',ScaledRawData,Model)%look at color overlay


%% comparison the results of "PR produced Zernike image" with "direct Zernike search"

clear all
clc

NA=1.49;lambda=0.69;n_med=1.33;n_imm=1.52;R=128; 

load('Y:\Farzin\Sequential Imaging Data\2015_7_30_CRLBpaper_test2_Collar19end\PR_result\darkredbead3_1-2015-7-30-16-52-27_PR_optim_single.mat')

OTFparam1=[1,obj.Results.SigmaX,obj.Results.SigmaY,0];
Z_start=obj.ParamsGui.Zrange.Value(1);
Z_end=obj.ParamsGui.Zrange.Value(3);
Zstack=[Z_start:0.1:Z_end];

dipshow(obj.Initial.ModifiedPsfData.plane1.psf)
dipshow(obj.Initial.PhaseRetrieveData.plane1.psf)
dipshow(obj.Initial.ZernikeFitData.plane1.psf) % Zernike produced image
dipshow(obj.Mpsf.plane1.Mpsf1)
PhaseCoeff=obj.Results.PhaseCoeff;
MagCoeff=obj.Results.MagCoeff;
NZ=7;%the spectrum of "n".                                       
%Theta0=[NA,lambda,n_med,n_imm,OTFparam1(1),PhaseCoeff,MagCoeff];
zernike_psf = Zernike_PSF(PhaseCoeff,MagCoeff,Zstack);
dipshow(zernike_psf)

joinchannels('RGB',zernike_psf,obj.Initial.ZernikeFitData.plane1.psf)



%% fitting to the OTF / write a class for finding NA from data
load('Y:\Farzin\Sequential Imaging Data\2015_7_6_CRLBbeads\PR_result\darkredbead-2015-7-6-16-5-46_PR_optim_single.mat')
a=obj.Data;
dipshow(obj.Data)
FOV=128;
PSF=obj.Data(:,:,11);% at the focal plane 
u1=fft2(PSF);
OTF=fftshift(u1); % look at the OTF in log stretch
%l=log(OTF);
p=OTF.^2;
u2=p(:,FOV/2-1);
%u2=OTF(:,FOV/2-1);
u3=abs(u2);
%u4=radialmean(u3) % doesn't look good 
u4=u3(65:127,:);
u5=sqrt(u4); % why sqrt: u4 changes fast and a lot, but sqrt(u4) will be less steep
lim1=1*max(u5);
lim2=0.05*max(u5);
mask=u5>lim2 & u5<lim1;
u6=u5(mask);
x=[1:FOV];
x1=x(mask);
x2=x1';
T=table(x2,u6);
figure
plot(x2,u6,'o')
[p,~,mu]=polyfit(T.x2,T.u6,3)
f=polyval(p,x2,[],mu)
hold on
plot(x2,f)
hold off
x3=[0:FOV/2-1];
f1=polyval(p,x3,[],mu)
hold on
plot(x3,f1)
[f1_min,x3_min]=min(f1)
FOV=128; % Field of view: total number of pixels in real space
bpp=96; % back projected pixel size (nm)
% lambda=690; % emission wavelength(nm)
delta_k=1/(FOV*bpp); % unit size in k-space
R=x3_min*delta_k; % radius of OTF (it's in k-space)
lambda_emm=690;
NA_over_lambda=R/2;
NA=NA_over_lambda*lambda_emm
% to show there is an error
delta1=0.75*max(f1)-lim1; 
delta2=0.25*max(f1)-lim2;

