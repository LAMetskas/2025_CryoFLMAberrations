function [psf_biplane] = genpsf_biplane_real(x,w,obj,dist_biplane)
R = obj.Boxsize;
obj.Xpos = x(:,1).*w(1);
obj.Ypos = x(:,2).*w(2);
Zpos_change = [-dist_biplane/2 dist_biplane/2];

N = numel(x(:,1));
% I = x(:,4);
% tmp = zeros(1,1,N);
% tmp(1,1,:) = I;
% IL = repmat(tmp,[R,R,1]);


psf_biplane = zeros(R,R,N,2);
for nn = 1:2
    obj.Zpos = x(:,3).*w(3) + Zpos_change(nn);

    obj.precomputeParam();    
    obj.padPupil();
    obj.genPSF_2();    
%     dipshow(obj.PSFs)
    
%     obj.genPupil();
%     obj.genPSF();
    obj.scalePSF('normal');
    psfI = obj.ScaledPSFs;
    
%     psfI = obj.PSFs; 
%     normf = sum(sum(obj.Pupil.mag.^2,1),2);
%     psf = psfI./normf;
%     psf = psfI./sum(psfI(:));
    normf = squeeze(sum(sum(psfI, 1), 2));  % size = [1󪻖01] ? [1�601]
    normf = reshape(normf, [1 1 N]);  % Nzs = 601
    psf = psfI ./ normf;

%     psf = psfI./normf;
    

    
    I = x(:,nn+3);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = I;
    IL = repmat(tmp,[R,R,1]).*w(4);

    bg = x(:,nn+5);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = bg;
    bgL = repmat(tmp,[R,R,1]).*w(5);
    
    psf_biplane(:,:,:,nn) = psf.*IL+bgL;
end

end
