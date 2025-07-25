function genZKresult(obj)
pupil_phase=obj.PRstruct.Pupil.phase;
pupil_mag=obj.PRstruct.Pupil.mag;
n=obj.RefractiveIndex;
Freq_max=obj.NA/obj.Lambda;
R=obj.PSFsize;
N=obj.DatadimZ;

NA_constrain=obj.k_r<Freq_max;
R_aber=angle(pupil_phase).*NA_constrain;
U=pupil_mag.*cos(R_aber).*NA_constrain;
V=pupil_mag.*sin(R_aber).*NA_constrain;
complex_Mag=U+1i*V;

[CN_complex,pupil_complexfit]=obj.fitzernike(complex_Mag,'mag');
[CN_phase,pupil_phasefit]=obj.fitzernike(R_aber,'phase');
[CN_mag,pupil_magfit]=obj.fitzernike(pupil_mag,'mag');

pupilcompare=pupil_complexfit-pupil_magfit.*exp(pupil_phasefit.*1i);
% generate Zernike fitted psf
zernike_psf=zeros(R,R,N);
z=[obj.Zstart:obj.Zstep:obj.Zend];
k_z=sqrt((n/obj.Lambda)^2-obj.k_r.^2).*NA_constrain;
for j=1:N
    defocus_phase=2*pi*z(j).*k_z.*1i;
    %pupil_complex=pupil_complexfit.*exp(defocus_phase);
    pupil_complex=pupil_magfit.*exp(pupil_phasefit.*1i).*exp(defocus_phase);
    psfA=abs(fftshift(fft2(pupil_complex)));
    Fig2=psfA.^2;
    zernike_psf(:,:,j)=Fig2./R^4;
end

obj.PSFstruct.ZKpsf=zernike_psf;
obj.PRstruct.Zernike_phase=CN_phase;
obj.PRstruct.Zernike_mag=CN_mag;
obj.PRstruct.Zernike_complex=CN_complex;
obj.PRstruct.Fittedpupil.complex=pupil_complexfit;
obj.PRstruct.Fittedpupil.phase=pupil_phasefit;
obj.PRstruct.Fittedpupil.mag=pupil_magfit;
end
