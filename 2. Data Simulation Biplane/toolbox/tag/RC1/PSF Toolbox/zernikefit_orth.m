
function [pupil_phaseZ,pupil_magZ,zernike_psf,residue_phase,residue_mag,Coeff_phase,Coeff_mag]=zernikefit_orth(pupilphase,pupilmag,Zstack,NZ,R,lambda,NA,pixelsizeCCD,Magnify,n_imm)

scale=R*pixelsizeCCD/Magnify;
k_r=double(rr(R,R)./scale);
pupil_phaseZ=zeros(R,R);
pupil_magZ=zeros(R,R);
n=n_imm;
Freq_max=NA/lambda;
NA_constrain=k_r<Freq_max; % we are making a mask 
R_aber=pupilphase.*NA_constrain; % left and right sides of equation are each 128-by-128 matrices
%NOTE: here we only use the in-focus pupilphase so it's not 128-128-Zstack but its only 128-by-128
Factor=zeros(1,2*NZ+2); % NZ is not the number of Zernike coefficients. It just shows "n"

for ii=1:2*NZ+2
    Factor(ii)=factorial(ii-1);
end

% generate Zernike terms in equation 56(128-by-128-by-36: all zernike coeffs for each pixel of the "in-focus" plane)
[Z]=zernike_terms1(NZ,R,lambda,NA,pixelsizeCCD,Magnify,Factor); % ? why just not Z instead of [Z] ?

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
Coeff_phase=CeffZ1; % a 1-by-36 matrix

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
    
    pupil_phaseZ=pupil_phaseZ+Z(:,:,k).*Coeff_phase(k);
    
end

residue_phase=R_aber-pupil_phaseZ;

for k=1:Zterms
    
    pupil_magZ=pupil_magZ+Z(:,:,k).*Coeff_mag(k);
    
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

