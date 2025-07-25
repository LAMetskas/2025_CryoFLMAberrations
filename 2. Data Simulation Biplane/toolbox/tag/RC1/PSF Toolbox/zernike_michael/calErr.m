function [sse,LL]=calErr(data1,data2,model1,model2)
R1=31;
szd=size(data1);
szm=size(model1);

realsize0=floor(R1/2);
realsize1=ceil(R1/2);
sz=szd;
starty=-realsize0+floor(sz(2)/2)+1;endy=realsize1+floor(sz(2)/2);
startx=-realsize0+floor(sz(1)/2)+1;endx=realsize1+floor(sz(1)/2);
data1i=data1(starty:endy,startx:endx,:);
data2i=data2(starty:endy,startx:endx,:);
sz=szm;
starty=-realsize0+sz(2)/2+1;endy=realsize1+sz(2)/2;
startx=-realsize0+sz(1)/2+1;endx=realsize1+sz(1)/2;
model1i=model1(starty:endy,startx:endx,:);
model2i=model2(starty:endy,startx:endx,:);
sz=size(model1i);
model1o=zeros(sz);
model2o=zeros(sz);

I=1000;
bg=2;
for ii=1:sz(3)
    est1=fminsearch(@(x)modelFit(x,model1i(:,:,ii),data1i(:,:,ii)),[I,bg],optimset('MaxIter',50,'Display','off'));
    est1(est1<0)=0;
    model1o(:,:,ii)=model1i(:,:,ii).*est1(1)+est1(2);
    est2=fminsearch(@(x)modelFit(x,model2i(:,:,ii),data2i(:,:,ii)),[I,bg],optimset('MaxIter',50,'Display','off'));
    est2(est2<0)=0;
    model2o(:,:,ii)=model2i(:,:,ii).*est2(1)+est2(2);
end

overlay=joinchannels('RGB',data1i,model1o);
sse=sum(sum(sum((data1i-model1o).^2+(data2i-model2o).^2)));
LL=sum(sum(sum(2*(data1i-model1o-model1o.*log(data1i)+model1o.*log(model1o)+data2i-model2o-model2o.*log(data2i)+model2o.*log(model2o)))));

end

function sse=modelFit(x,model,data)
I=x(1);
bg=x(2);
modeli=model.*I+bg;
sse=sum(sum((modeli-data).^2));
end

