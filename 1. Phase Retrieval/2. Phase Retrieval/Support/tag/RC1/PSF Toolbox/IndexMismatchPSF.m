function [psf,aber_phase]=IndexMismatchPSF(Refindex,NA,lambda,R,boxsize,pixelSize,lateralM,Zstack,depth,OTFparam,I,bg)
% Note: here we are simulating the beads on the coverglass from Sheng's paper "3D single molecule localization using PR function"
% but we might need to also have simulation for bead away from coverglass
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
deltaH=depth*n_med.*cos_theta2-depth*n_imm^2/n_med.*cos_theta1; % the last term of equ(9) is used at line 25
aber_index=2*pi/lambda.*deltaH.*NA_constrain; % aberration term
aber_phase=exp(aber_index.*1i).*NA_constrain;
% finding PSF by using the Pupil-Functuin of each plane of the z-stack:
for ii=1:length(Zstack)
    beta3=2*pi/lambda*n_imm*d3(ii).*cos_theta1.*NA_constrain; % defocus term
    defocus_phase=exp(1i.*beta3).*NA_constrain;
    pupil_phase=aber_phase.*defocus_phase; % pupil plane's phase. But later only the aber_phase is used to find Zernike Coeffs
    % The pupil_phase here is different for different z planes in z-stack.
    % but when calculating the zernike coeffs, we just use the focused plane
    PSF_A=ft(NA_constrain.*pupil_phase); % PSF_A=ft(NA_constrain*Pupil_Phase)
    psf0=abs(PSF_A).^2; % PSF  
    
    % OTF rescaling (we want the sharp edges to be smooth):
    OTF=OTFparam(1).*exp(-k_x.^2./(2*OTFparam(2)^2)).*exp(-k_y.^2./(2*OTFparam(3)^2))+OTFparam(4); % This is a smoothing filter in frequency space
    tp2=ift(psf0).*OTF; % we are miltiplying a filter named OTF by real otf which is ift(psf0) 
    tp3=abs(ft(tp2)); % after smoothing in frequency space, we go back :)
    Mod_psf=tp3./sum(tp3); % modified PSF
    start_y=-boxsize/2+R/2;end_y=boxsize/2+R/2-1; % 0 to 127
    start_x=-boxsize/2+R/2;end_x=boxsize/2+R/2-1; % 0 to 127
    model=I.*Mod_psf+bg; % this is our "model" : I*mod_PSF+bg
    Bipsf=model(start_y:end_y,start_x:end_x); % note that we use a single plane microscope not biplane
    psf(:,:,ii)=double(Bipsf);
end
end