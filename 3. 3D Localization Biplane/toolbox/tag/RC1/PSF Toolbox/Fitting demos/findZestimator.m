function obj=findZestimator(obj,data1,data2)
% findZestimator - find coefficients for initial estimation of z positions
% in 3D localization. 
%   The initial z estimation is found by the following steps:
%   1. estimate the total photon count of the measured PSFs from plane A and
%      plane B, denote as Ip1 and Ip2.
%   2. calculate Er = Ip2./Ip1.
%   3. find coefficients of linear fit of Er relative to z
%   4. save the coefficients in 'Zestimator' for 3D localization

zL=obj.Zpos_calib;
fitzRg=obj.FitZrange;
mask=(zL>=fitzRg(1))&(zL<=fitzRg(2));
data1i=data1(:,:,mask);
data2i=data2(:,:,mask);
boxsize=obj.BoxSizeFit;
Start=floor(boxsize/2);
s=obj.CenterSize;
z=zL(mask);
Num=length(z);
Ip1=zeros(Num,1);
Ip2=zeros(Num,1);
for ii=1:Num
    Row1=[mean(data1i(1,:,ii)),mean(data1i(boxsize,:,ii)),mean(data1i(:,1,ii)),mean(data1i(:,boxsize,ii))];
    Row2=[mean(data2i(1,:,ii)),mean(data2i(boxsize,:,ii)),mean(data2i(:,1,ii)),mean(data2i(:,boxsize,ii))];
    bg1=min(Row1);
    bg2=min(Row2);
    bg=min([bg1,bg2]);
    Ip1(ii)=sum(sum(data1i(Start+1-s:Start+1+s,Start+1-s:Start+1+s,ii)-bg));
    Ip2(ii)=sum(sum(data2i(Start+1-s:Start+1+s,Start+1-s:Start+1+s,ii)-bg));
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