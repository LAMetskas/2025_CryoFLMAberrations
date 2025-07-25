function obj=fitPSFdata_Bi(obj,data1,data2,coords1,coords2)
% fitPSFdata - find 3D localization results of the measured PSFs that were
% used for phase retrieval. 
%   This is used to test phase retrieved PSF model in 3D localization. It
%   uses mexfunction cudaBiplane_SM_Streams, which uses interpolation from
%   sample PSFs to generate the PSF model. Currently, it is only compatible
%   with EMCCD.
Num=size(data1,3);
maxFN=1e3;
m=ceil(Num/maxFN);
T=[reshape(repmat([1:m-1],maxFN,1),maxFN*(m-1),1);ones(Num-maxFN*(m-1),1).*m];
P=[];
LL=[];
SSE=[];
crlb=[];
CG=[];
% x=obj.PSFobj.Xpos;
% y=obj.PSFobj.Ypos;
% z=obj.PSFobj.Zpos;
%xIn=cat(2,x,y,Photon,Bg.*ones(Num,1),z);
for ii=1:m
    mask=(T==ii);
    %xIni=xIn(mask,:);
    data1i=squeeze(data1(:,:,mask));
    data2i=squeeze(data2(:,:,mask));
    coords1i=squeeze(coords1(mask,:));
    coords2i=squeeze(coords2(mask,:));
    [Pi,LLi,SSEi,crlbi,CGi]=PRfit(obj,data1i,data2i,coords1i,coords2i);
    P=cat(1,P,Pi);
    LL=cat(1,LL,LLi);
    SSE=cat(1,SSE,SSEi);
    crlb=cat(1,crlb,crlbi);
    CG=cat(1,CG,CGi);
end
obj.Fitresult.P=P;
obj.Fitresult.LL=LL;
obj.Fitresult.SSE=SSE;
obj.Fitresult.CRLB=crlb;
obj.Fitresult.CG=CG;
end

function [P,LL,SSE,crlb,CG]=PRfit(obj,data1,data2,coords1,coords2)
boxsize=obj.BoxSizeFit;
p=single(obj.Zestimator');
Num=size(data1,3);
BoxCenters1=cat(2,coords1,ones(Num,1));
BoxCenters2=cat(2,coords2,ones(Num,1));


channel1=reshape(single(data1),Num*boxsize^2,1);
channel2=reshape(single(data2),Num*boxsize^2,1);
Coords1=reshape(single(BoxCenters1(1:Num,:)'),Num*3,1);
Coords2=reshape(single(BoxCenters2(1:Num,:)'),Num*3,1);

[Pi,CGi,CRLBi,Erri,psf]=cudaBiplane_SM_Streams(channel1,channel2,...
    Coords1,Coords2,...
    obj.SamplePSF,p,obj.TransformP,...
    obj.SampleSpacingXY,obj.SampleSpacingZ,...
    obj.SampleS.StartX,obj.SampleS.StartY,obj.SampleS.StartZ,...
    obj.Iratio,obj.PlaneDis,obj.CenterSize,...
    obj.Iterationsfit,Num);
        

ResultSize=7;
x0Size=6;
P=reshape(Pi,ResultSize,Num)';
Err=reshape(Erri,2,Num)';
LL=Err(:,2);
SSE=Err(:,1);
crlb=reshape(CRLBi,x0Size,Num)';
CG=reshape(CGi,x0Size,Num)';  
PSF=reshape(psf,[16,16,2*Num]);
PSF1=PSF(:,:,1:2:2*Num-1);
PSF2=PSF(:,:,2:2:2*Num);
ov1=joinchannels('RGB',PSF1,data1);
ov2=joinchannels('RGB',PSF2,data2);

dipshow(ov1);
dipshow(ov2);

end





