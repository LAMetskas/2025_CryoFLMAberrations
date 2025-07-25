function [Zernc]=GenZCoeffv1(obj,sampleS,Thresh)

Nzern=size(obj.Z.ZM,3);
Z_N=sqrt(Nzern)-1;
pxsz=obj.Pixelsize;
psz=obj.PSFsize;


R_aber=obj.Pupil.phase;
U=obj.Pupil.mag.*cos(R_aber);
V=obj.Pupil.mag.*sin(R_aber);
complex_Mag1=U+1i*V;

CN_complex1=obj.Z.fitzernike(complex_Mag1,'mag',Z_N, psz);
g=pxsz*1e3/sampleS.Devx;
CN_complex=CN_complex1.*g;


[ZernC1,IndM1,IndN1,IndZern1]=ZernInd(Thresh,Nzern,CN_complex);
pCZ1_real=single(real(ZernC1)');
pCZ1_imag=single(imag(ZernC1)');
Zernc.pCZ1_real=pCZ1_real;
Zernc.pCZ1_imag=pCZ1_imag;
Zernc.IndM1=IndM1;
Zernc.IndN1=IndN1;
Zernc.IndZern1=IndZern1;

end
function [ZernC1,IndM1,IndN1,IndZern1]=ZernInd(Thresh,ObjectN,CeffCom1)
NZ=sqrt(ObjectN)-1;
Temp=abs(CeffCom1(1:ObjectN));
mask1=Temp>=Thresh;
ZernC1=CeffCom1(mask1);
j=2;
IndM=zeros(1,ObjectN);
IndN=zeros(1,ObjectN);
IndM(1)=0;
IndN(1)=0;
for nn=1:NZ
    for m=(j-1)/nn:-1:0
        if m~=0 
        IndM(j)=m;IndN(j)=nn; j=j+1;
        IndM(j)=m;IndN(j)=nn; j=j+1;
        else IndM(j)=m;IndN(j)=nn; j=j+1;
        end
    end
end 

IndZern=[0:1:ObjectN-1];
IndM1=single(IndM(mask1)');
IndN1=single(IndN(mask1)');
IndZern1=single(IndZern(mask1)');
end


