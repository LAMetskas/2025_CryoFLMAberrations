function fitPSFdata(obj)
% fitPSFdata - find 3D localization results of the measured PSFs that were
% used for phase retrieval. 
%   This is used to test phase retrieved PSF model in 3D localization. It
%   uses mexfunction cudaBiplane_SM_Streams, which uses interpolation from
%   sample PSFs to generate the PSF model. Currently, it is only compatible
%   with EMCCD.
tmp=load(fullfile(obj.PRobjA.FileDir,obj.PRobjA.FileName),'dataset');
namei=fields(tmp);
data=single(tmp.(namei{1}));
sz=size(data);
z=[obj.PRobjA.PRobj.Zstart:obj.PRobjA.PRobj.Zstep:obj.PRobjA.PRobj.Zend];
mask=(z>=obj.PRobjA.FitZrange(1))&(z<=obj.PRobjA.FitZrange(2));
ind=find(mask==1);
if numel(obj.PRobjA.PRobj.CCDoffset)>1
    CCDoffset=obj.PRobjA.PRobj.CCDoffset(:,:,1);
    Gain=obj.PRobjA.PRobj.Gain(:,:,1);
    CCDoffsetL=repmat(CCDoffset,[1,1,sz(3),sz(4)]);
    GainL=repmat(Gain,[1,1,sz(3),sz(4)]);
    data=(data-CCDoffsetL)./GainL;
else
    data=(data-obj.PRobjA.PRobj.CCDoffset)./obj.PRobjA.PRobj.Gain;
end
data=single(permute(data,[1,2,3,4]));

data(data<=0)=0.01;
P=[];
LL=[];
SSE=[];
crlb=[];
CG=[];
for ii=1:numel(ind)
    datai=squeeze(data(:,:,:,ind(ii)));
    [Pi,LLi,SSEi,crlbi,CGi]=PRfit(obj,datai);
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
function [P,LL,SSE,crlb,CG]=PRfit(obj,dataset)
boxsize=obj.BoxSizeFit;
data0=obj.PRobjA.PRobj.BeadData;
dimz0=obj.PRobjA.PRobj.DatadimZ;
p=single(obj.Zestimator');
sz=size(dataset);
if numel(sz)==3
    obj.PRobjA.PRobj.DatadimZ=sz(3);
    obj.PRobjB.PRobj.DatadimZ=sz(3);
else
    obj.PRobjA.PRobj.DatadimZ=1;
    obj.PRobjB.PRobj.DatadimZ=1;
end
[subregion1,Num]=genSubregion(obj.PRobjA,dataset);
[subregion2,Num]=genSubregion(obj.PRobjB,dataset);

est1=[256,1,1,1,0];
BoxCenters1=zeros(Num,3);
BoxCenters2=zeros(Num,3);
tmp1=repmat([80,60],Num,1);
[sse,tmp2]=TransformM1(est1,tmp1,tmp1);
BoxCenters1(:,1:2)=tmp1-floor(boxsize/2);
BoxCenters2(:,1:2)=tmp2-floor(boxsize/2);
TransformP=single(est1');

channel1=reshape(single(subregion1),Num*boxsize^2,1);
channel2=reshape(single(subregion2(:,:,1:Num)),Num*boxsize^2,1);
Coords1=reshape(single(BoxCenters1(1:Num,:)'),Num*3,1);
Coords2=reshape(single(BoxCenters2(1:Num,:)'),Num*3,1);

switch obj.FitType
    case 1
        [Pi,CGi,CRLBi,Erri,psf]=cudaBiplane_SM_Streams(channel1,channel2,...
                    Coords1,Coords2,...
                    obj.PRobjA.SamplePSF,p,TransformP,...
                    obj.PRobjA.SampleSpacingXY,obj.PRobjA.SampleSpacingZ,...
                    obj.PRobjA.SampleS.StartX,obj.PRobjA.SampleS.StartY,obj.PRobjA.SampleS.StartZ,...
                    obj.Iratio,obj.PlaneDis,obj.CenterSize,...
                    obj.Iterationsfit,Num);
        
%     case 2 % need to update the mexfunction
%         gainR=repmat(obj.GainRatio,[1,1,sz(3)]);
%         obj.PRobj.BeadData=gainR;
%         obj.PRobj.datapreprocess();
%         gainR1=obj.PRobj.Mpsf_subroi;
%         [gainR2]=chooseSubRegion(Dim1,Dim1,[1:1:Num],boxsize,gainR1);
%         gainRi=reshape(single(gainR2),Num*boxsize^2,1);
%         [Pi,CGi,CRLBi,Erri,psf]=cudaAst_SM_Streams(channel1,...
%             Coords1,...
%             obj.SamplePSF,...
%             obj.SampleSpacingXY,obj.SampleSpacingZ,...
%             obj.SampleS.StartX,obj.SampleS.StartY,obj.SampleS.StartZ,...
%             obj.Iterationsfit,Num,obj.FitType,x0i,gainRi);
end

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
maxf=max(max(max(mean(PSF1,3))),max(max(mean(subregion1,3))));
ov1=joinchannels('RGB',mean(PSF1,3)./maxf,mean(subregion1,3)./maxf);
ov1=single(ov1);
maxf=max(max(max(mean(PSF2,3))),max(max(mean(subregion2,3))));
ov2=joinchannels('RGB',mean(PSF2,3)./maxf,mean(subregion2,3)./maxf);
ov2=single(ov2);

figure('position',[400,100,600,400]);
subplot(231)
image(mean(subregion1,3),'cdatamapping','scale');
axis equal;axis tight;
subplot(232)
image(mean(PSF1,3),'cdatamapping','scale')
axis equal;axis tight;
subplot(233)
image(ov1);
axis equal;axis tight;
subplot(234)
image(mean(subregion2,3),'cdatamapping','scale');
axis equal;axis tight;
subplot(235)
image(mean(PSF2,3),'cdatamapping','scale')
axis equal;axis tight;
subplot(236)
image(ov2);
axis equal;axis tight;

colormap(jet)
% put back to origninal data
obj.PRobjA.PRobj.BeadData=data0;
obj.PRobjA.PRobj.DatadimZ=dimz0;
obj.PRobjB.PRobj.BeadData=data0;
obj.PRobjB.PRobj.DatadimZ=dimz0;

end

function [subregion,Num]=genSubregion(obj,dataset)
boxsize=obj.BoxSizeFit;
obj.PRobj.BeadData=dataset;
obj.PRobj.datapreprocess();
MpsfC1=obj.PRobj.Mpsf_subroi;
sz=size(MpsfC1);
Num=sz(3);
Dim1=ones(Num,1).*floor(sz(1)/2);
[subregion]=chooseSubRegion(Dim1,Dim1,[1:1:Num],boxsize,MpsfC1);
end



