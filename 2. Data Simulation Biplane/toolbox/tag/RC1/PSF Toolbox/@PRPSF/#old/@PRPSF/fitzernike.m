function [Ceffnorm,pupil_fitnorm]=fitzernike(obj,pupil,type)
R=obj.PSFsize;
Z=obj.Zk;
Ceff=[];
nz=0;
a=1;
b=1;
Zterms=(obj.ZernikeorderN+1)^2;
while nz<=obj.ZernikeorderN
    for ii=a:b
        temp=((2*nz+2)/pi)*sum(sum(pupil.*Z(:,:,ii)));
        Ceff(ii)=temp/((2*nz+2)/pi)/sum(sum(Z(:,:,ii).*Z(:,:,ii)));
    end
    nz=nz+1;
    a=a+2*nz-1;
    b=b+2*nz+1;
    if a>Zterms+1
        break; end
end

%generate fitted complex magnitude
pupil_fit=zeros(R,R);
for k=1:Zterms
    pupil_fit=pupil_fit+Z(:,:,k).*Ceff(k);
end

%normalize zernike coefficients
switch type
    case 'phase'
        tmp=exp(pupil_fit.*1i).*(1/R);
    case 'mag'
        tmp=pupil_fit.*(1/R);
end
tmp1=tmp.*conj(tmp);
tmp2=sum(sum(tmp1));
normF=sqrt(tmp2);
Ceffnorm=Ceff./normF;

pupil_fitnorm=zeros(R,R);
for k=1:Zterms
    pupil_fitnorm=pupil_fitnorm+Z(:,:,k).*Ceffnorm(k);
end
