function [z_guess,mcc_val,psf_model] = geniniBiplane_z_mat_parfor(data,psfobj_all,psftype,zLim,Nz,planedis)

R = size(data,1);
Nfit = size(data,3);
zi = linspace(zLim(1),zLim(2),Nz);
Nplane = size(data,4);
%psf_model = zeros(R,R,Nfit*Nz,Nplane);
psf_model = zeros(R,R,Nz,Nplane);  % simplified calculation, only cal once

for ss = 1:Nplane
    psfobj = psfobj_all{ss};
    psfobj.Xpos = zeros(Nz,1);% pixel
    psfobj.Ypos = zeros(Nz,1);% pixel
    switch psftype
        case 'normal'
            psfobj.Zpos = zi+planedis(ss);% micron
            psfobj.genPSF_2();
%             psfobj.genPSF();
            psfobj.scalePSF('normal');
%             psf_model(:,:,:,ss) = psfobj.PSFs;

        case 'IMM'
            psfobj.ZposMed = zi;% micron
            psfobj.genIMMPSF();
            psfobj.scalePSF('IMM');        
    end
    psf_model(:,:,:,ss) = psfobj.ScaledPSFs;
end
   

z_guess = zeros(Nfit,1);
mcc_val = zeros(Nfit,1);
img_fft2 = zeros(R,R,Nfit,Nplane);
ref_fft2 = zeros(R,R,Nz,Nplane);

for ii = 1: Nz
    %cal ref fft
    for ss = 1:Nplane
        ref = psf_model(:,:,ii,ss);
        ref = ref-mean(ref(:));
        ref = ref/std(ref(:));
        ref_fft2(:,:,ii,ss) = fft2(ref);
    end
end

parfor nn = 1:Nfit
   %cal img fft
   for ss = 1:Nplane
       img = data(:,:,nn,ss);
       img = img-mean(img(:));
       img = img/std(img(:));
       img_fft2(:,:,nn,ss) = fft2(img);
   end 
end

parfor nn = 1:Nfit
    mcc_max = 0;
    ind = 1;
    for ii = 1:Nz
        mcc = 0;
        for ss = 1:Nplane
            cc_value = abs(ifft2(ref_fft2(:,:,ii,ss) .*conj(img_fft2(:,:,nn,ss))));
            maxa = 1/(R*R)*max(cc_value(:));
            mcc = mcc + maxa/Nplane;
        end
        if mcc>mcc_max
            ind = ii;
            mcc_max = mcc;
        end
    end
    z_guess(nn) = zi(ind);
    mcc_val(nn) = mcc_max;
end

%%

% z_guess_2 = zeros(Nfit,1);
% mcc_val_2 = zeros(Nfit,1);
% 
% parfor nn = 1:Nfit
%     mcc_max = 0;
%     for ii = 1:Nz
%         mcc = 0;
%         for ss = 1:Nplane
%             ref = psf_model(:,:,(nn-1)*Nz+ii,ss);
%             img = data(:,:,nn,ss);
%             mcc = mcc + cc2(ref,img)/Nplane;
%         end
%         if mcc>mcc_max
%             ind = ii;
%             mcc_max = mcc;
%         end
%     end
%     z_guess_2(nn) = zi(ind);
%     mcc_val_2(nn) = mcc_max;
% end

end