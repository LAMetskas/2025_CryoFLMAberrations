function [Z]=zernike_terms1(NZ,R,lambda,NA,pixelsizeCCD,Magnify,Factor)
% let's work in polar coordinates (so we need r and phi for each point)
[X,Y]=meshgrid(-(-R/2:R/2-1),-R/2:R/2-1); % X & Y are both matrices of Y-by-X dimension
r1=sqrt(X.^2+Y.^2); % Y-by-X matrix
PHI=atan2(Y,X);  % for each point on the mesh 
scale=R*pixelsizeCCD/Magnify;
% Good Example to understand meshgrid: 
%xgv = [1 2 3];
%ygv = [1 2 3 4 5];
%[X,Y] = meshgrid(xgv, ygv)
Freq_max=NA/lambda;
k_r=r1./scale;
NA_constrain=k_r<Freq_max;
rou=k_r.*NA_constrain./Freq_max; % ? 
theta=PHI.*NA_constrain; % since all angles don't get into our objective

Zterms=(NZ+1)^2; %? NZ is giving us the spectrum of possible "n" / Zterms are the number of coeffs
Z=zeros(R,R,Zterms); % 128-128-36
Z(:,:,1)=1.*NA_constrain;
j=2;
for nn=1:NZ  % NZ=5
    for m=nn:-1:0
        RZ=newim(R,R);
        % calculating Q terms: (equation 55)
        for k=0:nn-m
            RZ=(-1)^k*Factor(2*nn-m-k+1)/(Factor(k+1)*Factor(nn-k+1)*Factor(nn-m-k+1)).*rou.^(2*(nn-m-k))+RZ;
            % this is a summation on RZ in each turn of the loop
           % NOTE: here the Factor(k+1)=k! (defined in zernikefit_orth.m) 
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
end
