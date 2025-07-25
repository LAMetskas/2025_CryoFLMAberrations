function fitFigure(obj)
% fitFigure - generate figures from localization results of measured PSF
% data. 
P=obj.Fitresult.P;
LL=obj.Fitresult.LL;

x=P(:,1);
y=P(:,2);
z=P(:,5);


pxsz=obj.PRobjA.PRobj.Pixelsize;
deviationx=(x-mean(x)).*pxsz.*1e3;
deviationy=(y-mean(y)).*pxsz.*1e3;

zrange=[obj.PRobjA.PRobj.Zstart,obj.PRobjA.PRobj.Zstep,obj.PRobjA.PRobj.Zend];
z1=[zrange(1):zrange(2):zrange(3)];
mask=(z1>=obj.PRobjA.FitZrange(1))&(z1<=obj.PRobjA.FitZrange(2));

zRef=z1(mask);
N=length(zRef);
Fn=length(x)/N;
vec=reshape(repmat([1:N],Fn,1),Fn*N,1);
deviationz=[];
for ii=1:N
    mask=vec==ii;
    deviationz=cat(1,deviationz,(z(mask)-zRef(ii)).*1e3);
end
zRefL=(vec-1)*zrange(2)+zrange(1);
ftsz=12;
figure('position',[100,200,800,700]);
subplot(311)
plot(zRefL,deviationx,'.');hold on;plot(zRef,mean(reshape(deviationx,Fn,N),1),'r-o')
xlabel('stage position (\mum)','fontsize',ftsz);ylabel('x deviaiton (nm)');
%ylim([-20,20]);xlim([zrange(1)-0.05,zrange(3)+0.05]);
subplot(312)
plot(zRefL,deviationy,'.');hold on;plot(zRef,mean(reshape(deviationy,Fn,N),1),'r-o')
xlabel('stage position (\mum)');ylabel('y deviaiton (nm)');
%ylim([-20,20]);xlim([zrange(1)-0.05,zrange(3)+0.05]);
subplot(313)
plot(zRefL,deviationz,'.');hold on;plot(zRef,mean(reshape(deviationz,Fn,N),1),'r-o')
xlabel('stage position (\mum)');ylabel('z deviaiton (nm)');
%ylim([-200,200]);xlim([zrange(1)-0.05,zrange(3)+0.05]);

% ImgName='XYZ-deviation-PRfit-ThreshZern-0d1-ROI50-k5-constantBg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-r300', '-djpeg', [obj.SaveDir,'\', ImgName])
%
boxsize=obj.PRobjA.BoxSizeFit;
X2_CDF=@(k,x)gammainc(x/2,k/2);
k=2*boxsize^2-6;
X2=LL;
pvalue=1-X2_CDF(k,X2);
figure('position',[100,200,800,700]);
subplot(211)
plot(zRef,mean(reshape(pvalue,numel(pvalue)/N,N),1),'r-o')
axis tight
xlabel('stage position (\mum)')
ylabel('pvalue')
subplot(212)
plot(zRef,mean(reshape(LL,numel(LL)/N,N),1),'r-o')
xlabel('stage position (\mum)')
ylabel('Log likelihood')
axis tight

end