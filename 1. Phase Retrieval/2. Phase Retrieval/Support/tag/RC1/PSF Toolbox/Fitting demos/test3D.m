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

data0=psfobj.ScaledPSFs;
Num=size(data0,3);
Bg=2;
mu=1000;
Photon=random('exp',mu,[Num,1]);
data=zeros(size(data0));
for ii=1:Num
data(:,:,ii)=data0(:,:,ii).*Photon(ii)+Bg;
end

data_noise=single(noise(data,'poisson',1));

obj.Iterationsfit=100;
obj=fitPSFdata_Ast(obj,data_noise,Bg,Photon);
%% plot fitting result
x0=obj.PSFobj.Xpos;
y0=obj.PSFobj.Ypos;
z0=obj.PSFobj.Zpos;

xfit=obj.Fitresult.P(:,1)-80;
yfit=obj.Fitresult.P(:,2)-60;
zfit=obj.Fitresult.P(:,5);

plot(z0-zfit)
hold on
plot(x0-xfit)
plot(y0-yfit)



