function obj=fitPSFdata_Ast(obj,data,Bg,Photon)
% fitPSFdata - find 3D localization results of the measured PSFs that were
% used for phase retrieval. 
%   This is used to test phase retrieved PSF model in 3D localization. It
%   uses mexfunction cudaAst_SM_Streams, which uses interpolation from sample
%   PSFs to generate the PSF model. It is compatible with both EMCCD and
%   sCMOS cameras.

Num=size(data,3);
maxFN=1e3;
m=ceil(Num/maxFN);
T=[reshape(repmat([1:m-1],maxFN,1),maxFN*(m-1),1);ones(Num-maxFN*(m-1),1).*m];
P=[];
LL=[];
SSE=[];
crlb=[];
CG=[];
x=obj.PSFobj.Xpos;
y=obj.PSFobj.Ypos;
z=obj.PSFobj.Zpos;
xIn=cat(2,x,y,Photon,Bg.*ones(Num,1),z);
for ii=1:m
    mask=(T==ii);
    xIni=xIn(mask,:);
    datai=squeeze(data(:,:,mask));
    [Pi,LLi,SSEi,crlbi,CGi]=PRfit(obj,datai,xIni);
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
function [P,LL,SSE,crlb,CG]=PRfit(obj,dataset,xIn)

boxsize=obj.BoxSizeFit;
Num=size(dataset,3);
BoxCenters1=zeros(Num,3);
tmp1=repmat([80,60],Num,1);
BoxCenters1(:,1:2)=tmp1-floor(boxsize/2)-0.5;

channel1=reshape(single(dataset),Num*boxsize^2,1);
Coords1=reshape(single(BoxCenters1(1:Num,:)'),Num*3,1);
%x0=genIniguess(dataset,z);
x0=xIn;
x0(:,1:2)=xIn(:,1:2)+floor(boxsize/2)+0.5;
x0i=single(reshape(x0',Num*5,1));

[Pi,CGi,CRLBi,Erri,psf]=cudaAst_SM_Streams(channel1,...
    Coords1,...
    obj.SamplePSF,...
    obj.SampleSpacingXY,obj.SampleSpacingZ,...
    obj.SampleS.StartX,obj.SampleS.StartY,obj.SampleS.StartZ,...
    obj.Iterationsfit,Num,obj.FitType,x0i);
        
 

ResultSize=6;
x0Size=5;
P=reshape(Pi,ResultSize,Num)';
Err=reshape(Erri,2,Num)';
LL=Err(:,2);
SSE=Err(:,1);
crlb=reshape(CRLBi,x0Size,Num)';
CG=reshape(CGi,x0Size,Num)';  
PSF=reshape(psf,[16,16,Num]);
ov=joinchannels('RGB',PSF,dataset);
dipshow(ov)
end


function x0=genIniguess(ROIStack1,z)
sz=size(ROIStack1);
Num=sz(3);
start=floor(sz(1)/2);
% xy
subroi=ROIStack1(start-2:start+2,start-2:start+2,:);
[x,y]=meshgrid([start-2:start+2],[start-2:start+2]);
xL=repmat(x,[1,1,Num]);
yL=repmat(y,[1,1,Num]);
area=squeeze(sum(sum(subroi,1),2));
comx=squeeze(sum(sum(subroi.*xL,1),2))./area;
comy=squeeze(sum(sum(subroi.*yL,1),2))./area;
%I,bg
subroi=ROIStack1(2:end-1,2:end-1,:);
sz1=size(subroi);
bg=(sum(sum(ROIStack1,1),2)-sum(sum(subroi,1),2))./(sz(1)^2-sz1(1)^2);
BG=repmat(bg,[sz(1),sz(2),1]);
I=squeeze(sum(sum(ROIStack1-BG,1),2));
I(I<200)=200;
bg=squeeze(bg);
x0=cat(2,comx,comy,I,bg,ones(Num,1).*z);
end

