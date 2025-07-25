classdef PRPSF < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Zstart = -1; %micron
        Zend = 1; %micron
        Zstep = 0.4; %micron
        Zindstart = 1; %index
        Zindend = 6;
        Zindstep = 1;
        Beadcenter;
        BeadXYshift;
        BeadData;
        DatadimX;
        DatadimY;
        DatadimZ;
        PSFsize;
        SubroiSize;
        Pixelsize; %micron
        IterationNum;
        IterationNumK;
        ZernikeorderN;
        OTFratioSize;
        NA;
        Lambda; %micrion
        RefractiveIndex;
        CCDoffset;
        Gain;
        FileDir;
        FileName;
    end
    
    properties (SetAccess = private, GetAccess = private)
        % precomputed images for k space operation
        PhiC; 
        ZoC;
        k_r;
        Phi;
        Zo;
        % precomputed factorials and zernike polynomials
        Factorial;
        Zk; 
    end
    % output parameters
    properties (SetAccess = private, GetAccess = public)
        Mpsf_subroi;
        Mpsf_extend;        
        PRstruct;
        PSFstruct;
    end
    methods
        function obj=PRPSF

        end
        function prepdata(obj)
            obj.PRstruct.NA=obj.NA;
            obj.PRstruct.Lambda=obj.Lambda;
            obj.PRstruct.RefractiveIndex=obj.RefractiveIndex;
            % input validation
            load(fullfile(obj.FileDir,obj.FileName),'dataset','Params');
            in=double(dataset);
            in=squeeze(mean(in,3));
            %CCDoffset=repmat(squeeze(mean(Params.Background,3)),[1,1,szB(3)]);
            %in=(in-CCDoffset).*obj.Gain;
            in=(in-obj.CCDoffset).*obj.Gain;
            f=min(in(:));
            in=in-f+1e-4;
            obj.BeadData=permute(in,[2,1,3]);
            [obj.DatadimX,obj.DatadimY,obj.DatadimZ]=size(in);
        end
        function precomputeParam(obj)
            [XC,YC]=meshgrid(-(-obj.DatadimX/2:obj.DatadimX/2-1),-obj.DatadimY/2:obj.DatadimY/2-1);
            obj.PhiC=atan2(YC,XC);
            obj.ZoC=sqrt(XC.^2+YC.^2);
            
            [X,Y]=meshgrid(-(-obj.PSFsize/2:obj.PSFsize/2-1),-obj.PSFsize/2:obj.PSFsize/2-1);
            obj.Zo=sqrt(X.^2+Y.^2);
            scale=obj.PSFsize*obj.Pixelsize;
            obj.k_r=obj.Zo./scale;
            obj.Phi=atan2(Y,X);
            
            f=zeros(1,2*obj.ZernikeorderN+2);
            for ii=1:2*obj.ZernikeorderN+2
                f(ii)=factorial(ii-1);
            end
            obj.Factorial=f;
            obj.Zk=obj.zernikepolynomial();
        end
        
        function datapreprocess(obj)
            realsize0=floor(obj.SubroiSize/2);
            realsize1=ceil(obj.SubroiSize/2);
            MpsfC=zeros(obj.SubroiSize,obj.SubroiSize,obj.DatadimZ);
            
            for ii=1:obj.DatadimZ
                tmp=fftshift(ifft2(squeeze(obj.BeadData(:,:,ii))));
                shiftphase=obj.ZoC./obj.DatadimX.*cos(obj.PhiC).*(obj.DatadimX/2-obj.Beadcenter(1)-1)-obj.ZoC./obj.DatadimY.*sin(obj.PhiC).*(obj.DatadimY/2-obj.Beadcenter(2)-1);
                tmp1=fft2(tmp.*exp(-2*pi.*shiftphase.*1i));
                startx=-realsize0+obj.DatadimX/2+1;endx=realsize1+obj.DatadimX/2;
                starty=-realsize0+obj.DatadimY/2+1;endy=realsize1+obj.DatadimY/2;
                tmp2=abs(tmp1(starty:endy,startx:endx));
                MpsfC(:,:,ii)=tmp2;
            end
            obj.Mpsf_subroi=MpsfC;
        end
        
        function genMpsf(obj)
            R1=obj.SubroiSize;
            R=obj.PSFsize;
            N=obj.DatadimZ;
            
            [X1,Y1]=meshgrid(-R1/2:R1/2-1,-R1/2:R1/2-1);
            circleF=sqrt(X1.^2+Y1.^2);
            circleF(circleF<=R1/2)=1;
            circleF(circleF>R1/2)=0;
            if R==R1
                circleF=ones(R,R);
            end
            realsize0=floor(R1/2);
            realsize1=ceil(R1/2);
            starty=-realsize0+R/2+1;endy=realsize1+R/2;
            startx=-realsize0+R/2+1;endx=realsize1+R/2;
            MpsfL=zeros(R,R,N);
            for ii=1:N
                Mpsfo=obj.Mpsf_subroi(:,:,ii);
                Edge=[mean(Mpsfo(1,:)),mean(Mpsfo(R1,:)),mean(Mpsfo(:,1)),mean(Mpsfo(:,R1))];
                bg=max(Edge);
                Fig2=(Mpsfo-bg);
                Fig2=Fig2.*circleF;
                Fig2(Fig2<=0)=0;
                tmp=zeros(R,R);
                tmp(startx:endx,starty:endy)=Fig2;
                Fig2=tmp;
                minimum0=0;
                Fig2(Fig2<=minimum0)=minimum0;
                MpsfL(:,:,ii)=Fig2./sum(sum(Fig2));
            end
            obj.Mpsf_extend=MpsfL;
        end
       
        function saveObj(obj,saveName)
            filename=obj.FileName;
            if isempty(obj.SaveDir)
                error('PRPSF:NoSaveDir','save directry is empty')
            else
                PRFileName=fullfile(obj.SaveDir,[filename saveName '.mat']);
                save(PRFileName,'obj');
            end
        end
        function calCRLB(obj)
        end
        
        function calibrategain(obj)
        end
    end
    
end

