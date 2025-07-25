function NA_data=findOTF_Radius(obj,Bpp,FOV,Lim1,Lim2)
%   finds the edge of Optical Transfer Function
%   This function takes in the infocus plane data of a bead 
%   alongside back projected pixel size and Field of View, and give back 
%   the radius of the OTF which is 2 times NA over wavelength of the
%   emission fluoroscent light

            PSF=obj.RawData; % we are looking at the focal plane, so average over all frames taken
            if ndims(PSF)==4 % ndim=length(size())
                fprintf('Error: image dimension is 4. It should be 3 or 2')
            end
            if ndims(PSF)==3
                u1=mean(PSF,3);
            end
            if ndims(PSF)==2
                u1=PSF;
            end
            u1=fft2(PSF); % could use "single" also and then "ft"
            OTF=fftshift(u1); % we need to shift to the center
            %OTF=ft(PSF);
            %u1=single(OTF);
            %dipshow(log(OTF))
            u2=u1(:,FOV/2-1);
            u3=radialmean(u2);
            dipshow(u3)
            u4=abs(u3);
            u5=sqrt(abs(radialmean(u2));
            dipshow(u5)
            Lim1=max(u5);
            Lim2=0.05*max(u5);
            mask=u5>Lim2 & u5<Lim1;
            u6=u5(mask);
            x=[1:FOV]; % Field of view = number of pixels on the image (128)
            x1=x(mask);
            x2=x1';
            T=table(x2,u6);
            [p,~,mu]=polyfit(T.x2,T.u6,3);
            f=polyval(p,x2,[],mu);
            x3=[0:50];
            f1=polyval(p,x3,[],mu);
            [f1_min,x3_min]=min(f1);
            delta_k=1/(FOV*Bpp); % unit size in k-space
            OTF_R=x3_min*delta_k; % radius of OTF (it's in k-space)/ OTF_R=2NA/lambda_emm
            NA_data=OTF_R*0.69/2
end

