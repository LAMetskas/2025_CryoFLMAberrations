function [sse,fitR]=TransformM1(x,position1,position2)
a=x(1);
b=x(2);
Mx=x(3);
My=x(4);
theta=x(5);

fitT(:,1)=Mx.*position1(:,1)+a;
fitT(:,2)=My.*position1(:,2)+b;

fitR(:,1)=fitT(:,1).*cos(theta)-fitT(:,2).*sin(theta);
fitR(:,2)=fitT(:,1).*sin(theta)+fitT(:,2).*cos(theta);
sse=sum(sum((fitR-position2).^2));
end
