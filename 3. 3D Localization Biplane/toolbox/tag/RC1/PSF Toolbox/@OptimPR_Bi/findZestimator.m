function findZestimator(obj)
% findZestimator - find coefficients for initial estimation of z positions
% in 3D localization. 
%   The initial z estimation is found by the following steps:
%   1. estimate the total photon count of the measured PSFs from plane A and
%      plane B, denote as Ip1 and Ip2.
%   2. calculate Er = Ip2./Ip1.
%   3. find coefficients of linear fit of Er relative to z
%   4. save the coefficients in 'Zestimator' for 3D localization

zL=[obj.PRobjA.PRobj.Zstart:obj.PRobjA.PRobj.Zstep:obj.PRobjA.PRobj.Zend];
fitzRg=obj.PRobjA.FitZrange;
mask=(zL>=fitzRg(1))&(zL<=fitzRg(2));
Mpsfo1=obj.PRobjA.PRobj.Mpsf_subroi(:,:,mask);
Mpsfo2=obj.PRobjB.PRobj.Mpsf_subroi(:,:,mask);
sz=size(Mpsfo1);
boxnum=sz(3);
Dim1=ones(boxnum,1).*floor(sz(1)/2);
boxsize=obj.BoxSizeFit;
[subregion1]=chooseSubRegion(Dim1,Dim1,[1:1:boxnum],boxsize,Mpsfo1);
[subregion2]=chooseSubRegion(Dim1,Dim1,[1:1:boxnum],boxsize,Mpsfo2);
Start=floor(boxsize/2);
s=obj.CenterSize;
z=zL(mask);
Er=[];
p=[];
Ip1=[];
Ip2=[];
for ii=1:length(z)
    Row1=[mean(subregion1(1,:,ii)),mean(subregion1(boxsize,:,ii)),mean(subregion1(:,1,ii)),mean(subregion1(:,boxsize,ii))];
    Row2=[mean(subregion2(1,:,ii)),mean(subregion2(boxsize,:,ii)),mean(subregion2(:,1,ii)),mean(subregion2(:,boxsize,ii))];
    bg1=min(Row1);
    bg2=min(Row2);
    bg=min([bg1,bg2]);
    Ip1(ii)=sum(sum(subregion1(Start+1-s:Start+1+s,Start+1-s:Start+1+s,ii)-bg));
    Ip2(ii)=sum(sum(subregion2(Start+1-s:Start+1+s,Start+1-s:Start+1+s,ii)-bg));
end
Er=Ip2./Ip1;
p = polyfit(Er,z,1);
obj.Zestimator=p;
figure('position',[200,200,400,300]);
f = polyval(p,Er);
plot(Er,z,'o',Er,f);
xlabel('I ratio');
ylabel('z position (\mum)')
legend('data','linear fit');
end