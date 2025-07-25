function genPRfigs(obj,ImgType)
% genPRfigs - generate figures of phase retrieval result, including various PSFs, pupil 
% functions, and plots of zernike coefficients. 
%
%   Input parameter: ImgType - type of image to be generated, select from
%   'PSF', 'pupil' and 'zernike'
% Assuming Mpsf is a 3D matrix with dimensions [rows, columns, slices]
[rows, cols, ~] = size(obj.Mpsf_extend);  % Get the size of the Mpsf matrix

% Calculate the radius as half of the number of rows or columns%
%radius = floor(min(rows, cols) / 2);  % Using floor to ensure an integer value

% Now use this dynamically calculated radius in your code
switch ImgType
    case 'PSF'
        z = obj.Zpos;
        zind = [obj.Zindstart:obj.Zindstep:obj.Zindend];
        L = length(zind);
        Mpsf = obj.Mpsf_extend;
        Modpsf = obj.PSFstruct.PRpsf;
        radius_measure = floor(size(obj.Mpsf_extend,1)/ 2);
        radius_fit=floor(size(Modpsf,1)/2);
        h1 = [];
        h4 = [];
        
        % positions [left, bottom, width, height]
        figure('Color',[1,1,1],'Name','measured and phase retrieved PSF at sampled z positions','Resize','on','Units','normalized','Position',[0.3,0.3,0.76,0.16])
        for ii = 1:L
            h1(ii) = subplot('position',[(ii-1)/(L+1),0.5,1/(L+1),1/2]);
            image(double(squeeze(Mpsf(:,:,zind(ii)))),'CDataMapping','scaled','Parent',h1(ii))
            text(3,3,[num2str(z(zind(ii)),3),'\mum'],'color',[1,1,1]);

            h4(ii) = subplot('position',[(ii-1)/(L+1),0,1/(L+1),1/2]);
            image(double(squeeze(Modpsf(:,:,zind(ii)))),'CDataMapping','scaled','Parent',h4(ii))
        end
        h1(ii+1) = subplot('position',[L/(L+1),0.5,1/(L+1),1/2]);
        image(double(permute(squeeze(Mpsf(radius_measure-10:radius_measure+10,radius_measure,:)),[2,1])),'CDataMapping','scaled','Parent',h1(ii+1))
        text(3,3,['x-z'],'color',[1,1,1]);
        h4(ii+1) = subplot('position',[L/(L+1),0,1/(L+1),1/2]);
        image(double(permute(squeeze(Modpsf(radius_fit-10:radius_fit+10,radius_fit,:)),[2,1])),'CDataMapping','scaled','Parent',h4(ii+1))
        colormap(jet)
        axis([h1,h4],'equal')
        axis([h1,h4],'off')

        saveas(gca,'measured and phase retrieved PSF.tif');


    % case 'pupil'
    %     % Define a larger figure size for better visibility
    %     figure('Color',[1,1,1],...
    %         'Name','Phase Retrieved and Zernike Fitted Pupil Function',...
    %         'Resize','on',...
    %         'Units','normalized',...
    %         'Position',[0.3,0.3,0.4,0.4]) % Adjusted figure size
    % 
    %     h1 = [];
    %     RC = 128;
    %     Rsub = 127;
    % 
    %     % Adjust the size and position of the first subplot
    %     h1(1) = subplot('Position',[0.05,0.1,0.4,0.8]); % Left subplot
    %     image(double(obj.PRstruct.Fittedpupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),...
    %         'CDataMapping','scaled',...
    %         'Parent',h1(1))
    %     text(3,8,'Zernike Pupil Mag','color',[1,1,1]);
    % 
    %     % Adjust the size and position of the second subplot
    %     h1(2) = subplot('Position',[0.55,0.1,0.4,0.8]); % Right subplot
    %     image(double(obj.PRstruct.Fittedpupil.phase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),...
    %         'CDataMapping','scaled',...
    %         'Parent',h1(2))
    %     text(3,8,'Zernike Pupil Phase','color',[0,0,0]);
    % 
    %     %colormap(gray)
    %     colormap("parula")
    %     axis(h1,'equal')
    %     axis(h1,'off')
    case 'pupil'
        % Define a larger figure size for better visibility
        figure('Color',[1,1,1],...
            'Name','Phase Retrieved and Zernike Fitted Pupil Function',...
            'Resize','on',...
            'Units','normalized',...
            'Position',[0.3,0.3,0.4,0.4]) % Adjusted figure size
    
        h1 = [];
        RC = 128;
        Rsub = 127;
    
        % Set desired color ranges
        mag_range = [0, 1];  % Example range for magnitude visualization
        phase_range = [-10, 10];  % Example range for phase visualization
    
        % Adjust the size and position of the first subplot
        h1(1) = subplot('Position',[0.05,0.1,0.4,0.8]); % Left subplot
        imagesc(double(obj.PRstruct.Fittedpupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)), mag_range);  % Set color range
        text(3,8,'Zernike Pupil Mag','color',[1,1,1]);
    
        % Adjust the size and position of the second subplot
        h1(2) = subplot('Position',[0.55,0.1,0.4,0.8]); % Right subplot
        imagesc(double(obj.PRstruct.Fittedpupil.phase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)), phase_range);  % Set color range
        text(3,8,'Zernike Pupil Phase','color',[0,0,0]);
    
        % Apply colormap and set properties
        colormap("parula")
        axis(h1,'equal')
        axis(h1,'off')
    
        % Add a horizontal color bar below the magnitude subplot
        colorbar(h1(1), 'southoutside', 'Ticks', mag_range, 'TickLabels', {'Min', 'Max'});
    
        % Add a vertical color bar for the phase subplot with range -20 to 20
        phase_cbar = colorbar(h1(2), 'eastoutside');  % Set color bar to the right side
        phase_cbar.Limits = phase_range;  % Set color limits for the color bar
        phase_cbar.Ticks = [-20, 0, 20];  % Define major tick marks
        phase_cbar.TickLabels = {'-20', '0', '20'};  % Custom tick labels
        phase_cbar.FontSize = 20;  % Increase font size for better visibility

    %case 'pupil'
        % pupil function at plane1
        %figure('Color',[1,1,1],'Name',' phase retrieved and Zernike fitted pupil function','Resize','on','Units','normalized','Position',[0.3,0.3,0.22,0.11])
        %h1=[];
        %RC=128;
        %Rsub=127;
        %h1(1)=subplot('Position',[0,0,1/2,1/2]);
        %image(double(obj.PRstruct.Fittedpupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(1))
        %text(3,8,['Zernike pupil mag'],'color',[1,1,1]);
        %h1(2)=subplot('Position',[0.51,0,1/2,1/2]);
        %image(double(obj.PRstruct.Fittedpupil.phase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(2))
        %text(3,8,['Zernike pupil phase'],'color',[0,0,0]);
        %colormap(gray)
        %axis(h1,'equal')
        %axis(h1,'off')
     
        saveas(gca,'phase retrieved and Zernike fitted pupil function.tif');   %edited by Fan Xu
    case 'zernike'
        PlotZernikeC(obj.PRstruct.Zernike_phase,'phase');
        saveas(gca,'phase_Zernike coefficinet.tif');   %edited by Fan Xu
        
        PlotZernikeC(obj.PRstruct.Zernike_mag,'magnitude');
        saveas(gca,'magnitude_Zernike coefficinet.tif');   %edited by Fan Xu
end
end

function PlotZernikeC(CN_phase,type)
nZ=length(CN_phase);
vec=linspace(max(CN_phase)-0.1,min(CN_phase)+0.1,8);
dinv=vec(1)-vec(2);
ftsz=12;
figure('position',[200,200,700,300])
plot(CN_phase,'o-')
text(nZ+5,vec(1),['x shift: ', num2str(CN_phase(2),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(2),['y shift: ', num2str(CN_phase(3),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(3),['z shift: ', num2str(CN_phase(4),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(4),['Astigmatism: ', num2str(CN_phase(5),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(5),['Astigmatism(45^o): ', num2str(CN_phase(6),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(6),['Coma(x): ', num2str(CN_phase(7),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(7),['Coma(y): ', num2str(CN_phase(8),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(8),['Spherical: ', num2str(CN_phase(9),'%.2f')],'fontsize',ftsz);
xlim([0,nZ+30])
ylim([min(CN_phase)-dinv/2,max(CN_phase)+dinv/2])
set(gca,'fontsize',ftsz)
xlabel('Zernike coefficient number','fontsize',ftsz)
ylabel('Value','fontsize',ftsz)
title(type)
% ImgName='_ZernCoeff';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-r300', '-djpeg', [PRtest.FileDir,PRtest.FileName(1:end-4),ImgName])

end
