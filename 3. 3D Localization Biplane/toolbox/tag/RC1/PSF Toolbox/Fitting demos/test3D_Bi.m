% find z estimator
Num=21;
planeDis=0.4;%micron
Zpos=linspace(-1,1,Num)'-planeDis/2;
psfobj.Xpos=zeros(Num,1);
psfobj.Ypos=zeros(Num,1);

psfobj.Zpos=Zpos;
psfobj.genPSF();
psfobj.scalePSF();
psf1=psfobj.ScaledPSFs;

psfobj.Zpos=Zpos+planeDis;
psfobj.genPSF();
psfobj.scalePSF();
psf2=psfobj.ScaledPSFs;

Bg=2;
Photon=1000;
data1=psf1.*Photon+Bg;
data2=psf2.*Photon+Bg;

obj.Zpos_calib=Zpos;
obj.FitZrange=[-0.8,0.7];% micron
obj.BoxSizeFit=16;
obj.CenterSize=4;
obj.PlaneDis=planeDis;
obj=findZestimator(obj,data1,data2);

%% localization
Num=100;
x0=80.*rand(Num,1);
y0=60.*rand(Num,1);
Zpos=linspace(-1,1,Num)';
obj.TransformP=single([64,1,1,1,0.01]');

tmp1=cat(2,x0,y0);
[sse,tmp2]=TransformM(obj.TransformP,tmp1,tmp1);
center1=tmp1-round(tmp1);
center2=tmp2-round(tmp2);
x1=center1(:,1);
y1=center1(:,2);
x2=center2(:,1);
y2=center2(:,2);

Coords1=tmp1-floor(obj.BoxSizeFit/2)-center1;
Coords2=tmp2-floor(obj.BoxSizeFit/2)-center2;


psfobj.Xpos=x1;
psfobj.Ypos=y1;
psfobj.Zpos=Zpos; 
psfobj.genPSF();
psfobj.scalePSF();
psf1=psfobj.ScaledPSFs;

psfobj.Xpos=x2;
psfobj.Ypos=y2;
psfobj.Zpos=Zpos+obj.PlaneDis; 
psfobj.genPSF();
psfobj.scalePSF();
psf2=psfobj.ScaledPSFs;

% dipshow(psf1)
% dipshow(psf2)

Bg=2;
Photon=1000;
data1=psf1.*Photon+Bg;
data2=psf2.*Photon+Bg;

%%
obj.BoxSizeFit=16;
obj.SampleS.Devz=0.005; %um
obj.SampleS.Dlimz=single([-1.5,1.5])'; % um
obj.SampleS.x0Size=5;
obj.SampleS.Zstack=[obj.SampleS.Dlimz(1):obj.SampleS.Devz:obj.SampleS.Dlimz(2)]';
obj.PSFobj=psfobj;
obj.OTFobj=OTFrescale();
obj.SampleS.PixelSize=obj.PSFobj.Pixelsize;
obj.SampleS.Devx=obj.SampleS.PixelSize/4*1e3; %nm
obj.SampleS.PsfSizeFine=round(2*obj.BoxSizeFit*obj.SampleS.PixelSize/obj.SampleS.Devx*1e3);
obj.SampleS.PixelSizeFine=2*obj.BoxSizeFit*obj.SampleS.PixelSize/obj.SampleS.PsfSizeFine;
w=obj.SampleS.PixelSize./obj.SampleS.PixelSizeFine;
sigma=[2,2].*w^2;
obj.SampleS.OTFparam1=single([1/w,sigma(1),sigma(2),0])';

obj.FitType=1;% for EMCCD
obj=genSamplePSF(obj);
obj.Iratio=1;
obj.Iterationsfit=100;
obj=fitPSFdata_Bi(obj,data1,data2,Coords1,Coords2);
%% plot fitting result
x0=x1;
y0=y1;
z0=Zpos;

xfit=obj.Fitresult.P(:,1)-0.5-round(obj.Fitresult.P(:,1)-0.5);
yfit=obj.Fitresult.P(:,2)-0.5-round(obj.Fitresult.P(:,2)-0.5);
zfit=obj.Fitresult.P(:,5);

plot(z0-zfit)
hold on
plot(x0-xfit)
plot(y0-yfit)

