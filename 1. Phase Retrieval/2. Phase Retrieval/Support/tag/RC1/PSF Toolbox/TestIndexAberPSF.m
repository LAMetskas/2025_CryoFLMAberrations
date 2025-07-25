%% generated PSF and phase for index mismatch aberration
addpath('C:\Users\Farzin\Dropbox\SIM\Isotropic CRLB\Sheng Codes')
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
OTFparam1=[1,2,2,0]; % [magnitude,sigmaX,sigmaY,bg]
R=128;
% IMPORTANT! depth: for beads stuck to coverslip should be zero. Note: if we have
% depth~=0, then if our model is correct we wouldn't have symmetry above and under zero in our synthetic PSF so it can work as a test of the model.
depth=single(0.0001); % micrion  
Zstack=[-1:0.1:1]';
% Important Note: here we are producing a PSF for known z-planes, but when
% we have the real data, we don't know what is the z of our bead exactly,
% so we don't have our pupil_phases at start and we need to find them in
% reverse and that's the PhaseRetrieval process (so this code is not the Phase Retrieval code!)

% here we find PSF by finding the Pupil-Function which itself is found by
% calculating two phase terms: 1- index mismatch phase  2- defocus phase
[PSF,aber_phase]=IndexMismatchPSF(Refindex,NA,lambda,R,boxsize,pixelsizeCCD,...
                    Magnify,Zstack,depth,OTFparam1,Intensity,BG1);

dipshow(PSF)
% angle(z)=atan2(imag(z),real(z))
pupilphase=angle(aber_phase); % "angle" returns phase angles for each element of complex array 
pupilmag=abs(aber_phase);
dipshow(pupilphase) % when i put depth=0, it gives a gray image, so i put depth=0.001.

%% expand in zernike polynomial
NZ=5;% Not the number of zernike coeffs but the spectrum of "n".
% we don't use the defocus phase in finding Zernike coeffs:
[pupil_phaseZ,pupil_magZ,zernike_psf,...
 residue_phase,residue_mag,Coeff_phase,Coeff_mag]=zernikefit_orth(pupilphase,pupilmag,...
                                            Zstack,NZ,R,lambda,NA,pixelsizeCCD,...
                                            Magnify,n_imm);
                                        
figure;plot(Coeff_phase)
title('Phase of Zernike Coefficients, depth=0')
xlabel('36 zernike terms')
ylabel('phase coefficient')
figure;plot(Coeff_mag)
title('Amplitude of Zernike Coefficients, depth=0')
xlabel('36 zernike terms')
ylabel('amplitude coefficient')

dipshow(pupil_phaseZ)
dipshow(zernike_psf)  % the only difference of zernike_psf with the PSF we calculated in the section above
% is that in zernike_psf we haven't used OTF-rescaling. If we do, the
% results would be the same (for NZ=10 maybe). Also when we use this
% zernike_psf, we can talk about our aberrations using zernike coeffs,
% which is a nice feature.

%% Minimizing the Mean Square Error:
% here we are minimizing the error between our Zernike_PSF model(above) and Raw Data 
% Goal: write a function as above, to do the MMSE.
% first load the PhaseRetrieved data of your bead:
addpath('Y:\Farzin\Sequential Imaging Data\2015_7_6_CRLBbeads\PR_result')
load('darkredbead-2015-7-6-16-5-46_PR_optim_single.mat')
RawData=obj.Data; % 128-128-21
error=zeros(R,R,length(Zstack));

for ii=1:length(Zstack)
% Note: zernike_psf is a psf produced synthetically
% Note: we need to rescale our RawData since it's intensity is 1000 times of zernike_psf data 
ScaledRawData(:,:,ii)=RawData(:,:,ii)/sum(sum(RawData(:,:,ii)));
Scaledzernike_psf(:,:,ii)=zernike_psf(:,:,ii)/sum(sum(zernike_psf(:,:,ii)));
error(:,:,ii)=Scaledzernike_psf(:,:,ii)-ScaledRawData(:,:,ii); % "error = model - data" at each pixel in each z-plane
S_error(:,:,ii)=error(:,:,ii).*error(:,:,ii)'; % square of error / 128-128-21
MSE=mean(mean(S_error));
% to find the MMSE, we need to calculate the covariance matrix of raw data:cov=E[(X-E(X))(X-E(X))]
end

model_MMSE=mean2(ScaledRawData)
dipshow(ScaledRawData)
dipshow(Scaledzernike_psf)
MSE=mse(Scaledzernike_psf,ScaledRawData) % This is just a scalar
dipshow(S_error)
dipshow(MSE)

% finding MMSE using fminsearch:
x1=ScaledRawData;
x2=Scaledzernike_psf;
x=[x1,x2];
MSE=@(x)mse(x(1),x(2));
x20(1)=0.1;x20(2)=0.1;x20(3)=0.1; % initial guess
p0=x20;
[x,fval,exitflag,output]=fminsearch(MSE,p0)

% general example of how to rescale RawData:
% Goal: we want for each frame, the sum over all pixel values =1
A=rand(2);
B=sum(sum(A));
ScaledA=A/B;
sum(sum(ScaledA)); % this is 1. 

%% fminsearch method

% practice test1: one input and one output
x=[-1:0.01:3];
x0=1;
fhandle=@function1;
[x,fmin,exitflag,output]=fminsearch(fhandle,x0)
disp('Minumun of the function is')
disp(fmin)
disp('... and is located at x = ')
disp(x)
disp('')

% practice test2: Rosenbrock banana 
x1=[-2:0.01:2];
x2=[-2:0.01:2]; 
x3=[-2:0.01:2];
x=[x1,x2,x3];
a=3;
banana=@(x)100*(x(2)-x(1).^2).^2+(1-x(1)).^2+a*sin(x(3));
x10=1.4;
x20=1.2;
x30=0.2;
p0=[x10,x20,x30];
[x,fval,exitflag,output]=fminsearch(banana,p0) % you can't put x10 and x20 intread of p0











