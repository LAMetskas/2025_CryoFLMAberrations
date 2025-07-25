function [Z]=zernikepolynomial(obj)
R=obj.PSFsize;
PHI=obj.Phi;
k_r=obj.k_r;
Factor=obj.Factorial;
Freq_max=obj.NA/obj.Lambda;
NA_constrain=k_r<Freq_max;
rou=k_r.*NA_constrain./Freq_max;
theta=PHI.*NA_constrain;

j=2;
Zterms=(obj.ZernikeorderN+1)^2;
Z=zeros(R,R,Zterms);
Z(:,:,1)=1.*NA_constrain;
for nn=1:obj.ZernikeorderN
    for m=(j-1)/nn:-1:0
        RZ=newim(R,R);
        for k=0:nn-m
            RZ=(-1)^k*Factor(2*nn-m-k+1)/(Factor(k+1)*Factor(nn-k+1)*Factor(nn-m-k+1)).*rou.^(2*(nn-m-k))+RZ;
        end
        
        if m~=0 
        Z(:,:,j)=RZ.*rou.^m.*cos(m.*theta); j=j+1;
        Z(:,:,j)=RZ.*rou.^m.*sin(m.*theta); j=j+1;
        else Z(:,:,j)=RZ.*rou.^m.*NA_constrain; j=j+1;
        end
    end
end 
end