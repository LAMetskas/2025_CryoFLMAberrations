classdef OTF_BandLimit < handle
    %OTF_BandLimit Finds the edge of Optical Transfer Function
    %   This function takes in the infocus plane data of a bead
    %   alongside back projected pixel size and Field of View, and give back
    %   the radius of the OTF which is 2 times NA over wavelength of the
    %   emission fluoroscent light
    
    properties
        RawData;       % Raw bead data
        PSF;
        ShowSteps = 1;
        Bpp;  %Back-projected pixel size (micron)
        OTF_Radius;
        FOV;  %Field Of View
        Lambda_emm = .69;   %micron
        OTF_Limit   %2*NA/lambda in micron^-1
        NA
        ShowFit =1;
        MinFit = .2;
        MaxFit = .5;
    end
    
    
    methods
        
        function NA_data=findOTF_Radius(obj)
            %   finds the edge of Optical Transfer Function
            
            switch ndims(obj.RawData) % better to use "switch" than "if"
                case 2 %Do nothing, already 2D
                    obj.PSF=single(obj.RawData);
                case 3 %Average 3rd dimension
                    obj.PSF=mean(single(obj.RawData),3);
                otherwise
                    error('Data must be 2D or 3D')
            end
            
            obj.PSF=single(backgroundoffset(obj.PSF));
            
            OTF=fftshift(fft2(obj.PSF)); % we need to shift to the center
            OTF_RadMean = radialmean(abs(OTF)); % ? why the size is 91? why not 64?

            OffSet = 3; %Ignore bump at zero from PSF background;
            LogRadMean = log(OTF_RadMean(OffSet:floor(min(size(OTF))/2)-1));
            dipshow(LogRadMean)
            
            LFilter = laplace(LogRadMean,3); %second derivative filter
            OTF_Limit_Pixels=max(find(maxima(LFilter)))+OffSet-1;
            
            Delta_k=1/(size(obj.PSF,1)*obj.Bpp); % unit size in k-space
            obj.OTF_Limit=OTF_Limit_Pixels*Delta_k; %(=2NA/lambda_emmision)
            
            obj.NA = obj.Lambda_emm*obj.OTF_Limit/2;
            
            figure
            line('XData', [OTF_Limit_Pixels OTF_Limit_Pixels], 'YData', [0 max(LogRadMean)], 'LineStyle', '--', ...
                'LineWidth', 1, 'Color','b');
            legend('OTF CutOff')
            
            circle1=rr<OTF_Limit_Pixels & rr>OTF_Limit_Pixels-1; % it's 256by256
            circle2=cut(circle1,[length(obj.PSF),length(obj.PSF)]); % resize to 128by128
            OTF_Sim_vs_Data=overlay(abs(log(OTF)).*20,circle2)
            

        end
        
    end
    
end

