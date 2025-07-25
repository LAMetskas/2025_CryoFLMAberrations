function [ zernike_psf ] = Zernike_PSF1(Theta0,Zstack) %Theta0 is the initial guess. We get it from PR results

M=[1,1];
Magnify=250; 
pixelsizeCCD=24;% micrion
boxsize=128;
n_med=1.33;
n_imm=1.52;
Intensity=2000;
BG1=0;
Refindex=[n_imm,n_med];
lambda=0.637;% micrion
NA=1.49;
OTFparam1=[1,2,2,0]; % [magnitude,sigmaX,sigmaY,bg]
R=128;
%Zstack=[-1:0.1:1]';
% let's work in polar coordinates (so we need r and phi for each point)
[X,Y]=meshgrid(-(-R/2:R/2-1),-R/2:R/2-1); % X & Y are both matrices of Y-by-X dimension
r1=sqrt(X.^2+Y.^2); % Y-by-X matrix
PHI=atan2(Y,X);  % for each point on the mesh 
scale=R*pixelsizeCCD/Magnify;
Freq_max=NA/lambda;
k_r=r1./scale;
NA_constrain=k_r<Freq_max; % we are making a mask 
rou=k_r.*NA_constrain./Freq_max; % ? 
theta=PHI.*NA_constrain; % since all angles don't get into our objective
NZ=5;
Zterms=(NZ+1)^2; %? NZ is giving us the spectrum of possible "n" / Zterms are the number of coeffs
Z=zeros(R,R,Zterms); % 128-128-36
Z(:,:,1)=1.*NA_constrain;
scale=R*pixelsizeCCD/Magnify;
k_r=double(rr(R,R)./scale);
pupil_phaseZ=zeros(R,R);
pupil_magZ=zeros(R,R);
n=n_imm;
lateralM=250;
pixelSize=24;
scale=R*pixelSize/lateralM; % R=128
Freq_max=NA/lambda;
psf=zeros([boxsize,boxsize,length(Zstack)]); % boxsize < R
d3=Zstack;
n_imm=Refindex(1);  % immersion medium (Oil)
n_med=Refindex(2); % in water (Cell)
k_x=abs(xx(R,R)./scale); % the reason we define k_x like this, is for FT
k_y=abs(yy(R,R)./scale); 
k_r=double(rr(R,R)./scale); % k_r is radial component of vector k
NA_constrain=double(k_r<Freq_max); % this is a mask 
sin_theta1=k_r.*lambda./n_imm; % NOTE: k is defined as 1/lambda_n
sin_theta2=n_imm./n_med.*sin_theta1; % Snell's law
cos_theta2=sqrt(1-sin_theta2.^2);
cos_theta1=sqrt(1-sin_theta1.^2);
% now we write the phases due to equation 9 of the paper:
depth=single(0.0001); % micrion
deltaH=depth*n_med.*cos_theta2-depth*n_imm^2/n_med.*cos_theta1; % the last term of equ(9) is used at line 25
% IMPORTANT! depth: for beads stuck to coverslip should be zero. Note: if we have
% depth~=0, then if our model is correct we wouldn't have symmetry above and under zero in our synthetic PSF so it can work as a test of the model.
aber_index=2*pi/lambda.*deltaH.*NA_constrain; % aberration term
aber_phase=exp(aber_index.*1i).*NA_constrain;
for ii=1:length(Zstack)
    beta3=2*pi/lambda*n_imm*d3(ii).*cos_theta1.*NA_constrain; % defocus term
    defocus_phase=exp(1i.*beta3).*NA_constrain;
    pupil_phase=aber_phase.*defocus_phase; % pupil plane's phase. But later only the aber_phase is used to find Zernike Coeffs
    % The pupil_phase here is different for different z planes in z-stack,
    % but when calculating the zernike coeffs, we just use the focused plane
end
% angle(z)=atan2(imag(z),real(z))
pupilphase=angle(aber_phase); % "angle" returns phase angles for each element of complex array 
pupilmag=abs(aber_phase);
R_aber=pupilphase.*NA_constrain; % left and right sides of equation are each 128-by-128 matrices
%NOTE: here we only use the in-focus pupilphase so it's not 128-128-Zstack but its only 128-by-128
Factor=zeros(1,2*NZ+2); % NZ is not the number of Zernike coefficients. It just shows "n"
for ii=1:2*NZ+2
    Factor(ii)=factorial(ii-1);
end
j=2;
for nn=1:NZ  % NZ=5
    for m=nn:-1:0
        RZ=newim(R,R);
        % calculating Q terms: (equation 55)
        for k=0:nn-m
            RZ=(-1)^k*Factor(2*nn-m-k+1)/(Factor(k+1)*Factor(nn-k+1)*Factor(nn-m-k+1)).*rou.^(2*(nn-m-k))+RZ;
            % this is a summation on RZ in each turn of the loop
           % NOTE: here the Factor(k+1)=k! (defined above) 
        end
        % calculating equ 56: (here we don't calculate the coefficients of each zernike term. We only calculate the terms.)
        if m~=0 
        % Note: odd terms have "cos" and even terms have "sin" so with
        % j=j+1 we are going from an odd term to an even term and vice versa.
        Z(:,:,j)=RZ.*rou.^m.*cos(m.*theta); j=j+1; 
        Z(:,:,j)=RZ.*rou.^m.*sin(m.*theta); j=j+1;
        else Z(:,:,j)=RZ.*rou.^m.*NA_constrain; j=j+1;
        end
    end
end 

% generate coefficient for phase (for each of 36 terms we find proper coeff_phase)
CeffZ1=[];
nz=0;
a=1; % this is the starting "a" in the loop
b=1; % as above
Zterms=(NZ+1)^2; % We have 36 Zernike terms (NZ=5)
while nz<=NZ
    for ii=a:b
        % good example of sum(sum(A)) would be for A=[1 2 3;4 5 6]
        temp=((2*nz+2)/pi)*sum(sum(R_aber.*Z(:,:,ii))); % ? Z(:,:,ii) shows zernike term #ii for all the 128-by-128 pixels of the in-focus plane. 
        % Q)shouldn't R_aber be also a function of z? A) No! here we only need the term with defocus=0.
        CeffZ1(ii)=temp/((2*nz+2)/pi)/sum(sum(Z(:,:,ii).*Z(:,:,ii)));
        % we expand the pupil phase function on the zernike basis
        % CeffZ is 1-by-36: where most of them are zero 
    end
    nz=nz+1;
    a=a+2*nz-1;
    b=b+2*nz+1;
    if a>Zterms+1
        break; end
end
 
Coeff_phase1=CeffZ1; % a 1-by-36 matrix
Coeff_phase2=Theta0(6:end); % a 1-by-36 matrix / initial 

% generate coefficient for pupil magnitude
CeffZ2=[];
nz=0;
a=1;
b=1;
R_mag=pupilmag.*NA_constrain;
while nz<=NZ
    for ii=a:b
        temp=(2*nz+2)/pi*sum(sum(R_mag.*Z(:,:,ii)));
        CeffZ2(ii)=temp/((2*nz+2)/pi)/sum(sum(Z(:,:,ii).*Z(:,:,ii)));
        % we expand the pupil magnitude function on the zernike basis
    end
    nz=nz+1;
    a=a+2*nz-1;
    b=b+2*nz+1;
    if a>Zterms+1
        break; end
end

Coeff_mag=CeffZ2; % a 1-by-36 matrix

%generate fitted phase and magnitude
for k=1:Zterms

    pupil_phaseZ=pupil_phaseZ+Z(:,:,k).*Coeff_phase2(k); % You should also be able to
    % produce results with Coeff_phase1(k) if you don't have any initial
    % guess for Theta0, so add that option.
end

residue_phase=R_aber-pupil_phaseZ;

for k=1:Zterms
    
    pupil_magZ=pupil_magZ+Z(:,:,k).*Coeff_mag2(k);
    
end

tmp=pupil_magZ.*exp(pupil_phaseZ.*1i).*(1/R);
tmp1=tmp.*conj(tmp);
tmp2=sum(sum(tmp1));
normF=sqrt(tmp2);
Coeff_mag=Coeff_mag./normF;
pupil_magZ=zeros(R,R);

for k=1:Zterms
    
    pupil_magZ=pupil_magZ+Z(:,:,k).*Coeff_mag(k);
    
end
residue_mag=pupilmag-pupil_magZ;

% generate Zernike fitted data
N=length(Zstack);
zernike_psf=zeros(R,R,N);
k_z=sqrt((n/lambda)^2-k_r.^2).*NA_constrain;
for j=1:N
    defocus_phase=2*pi*Zstack(j).*k_z.*1i;
    pupil_complex=pupil_magZ.*exp(pupil_phaseZ.*1i).*exp(defocus_phase);
    psfA=abs(fftshift(fft2(pupil_complex))); % fftshift: shifts zero-freq component to center of spectrum
    % psf_A=2D_FT(pupil_complex)
    Fig2=psfA.^2;
    zernike_psf(:,:,j)=Fig2./R^4; % 128-128-21
    % NOTE: Zernike_psf is an approximation of Phase_Retrieved PSF that we find from the data
    % the reason is we have noise so the coefficients are not guessed completely correct
end

%dipshow(zernike_psf)
end

