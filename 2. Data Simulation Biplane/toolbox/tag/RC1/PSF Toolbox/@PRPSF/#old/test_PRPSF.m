% test PRPSF

probj=PRPSF;
probj.CCDoffset=100;
probj.Gain=0.5;
probj.NA=1.49;
probj.Lambda=0.67;
probj.RefractiveIndex=1.52;
probj.Pixelsize=0.1;
probj.PSFsize=128;
probj.SubroiSize=90;
probj.OTFratioSize=60;
probj.ZernikeorderN=7;
probj.Zstart = -1.4; %micron
probj.Zend = 1.4; %micron
probj.Zstep = 0.4; %micron
probj.Zindstart = 1; %index
probj.Zindend = 8;
probj.Zindstep = 1;
probj.IterationNum=25;
probj.IterationNumK=5;

[FileName FileDir filterInd] = uigetfile('Y:\*.mat','Select file:','MultiSelect','off');
probj.FileDir=FileDir;
probj.FileName=FileName;
probj.prepdata();
probj.precomputeParam();
probj.findXYshift();
probj.datapreprocess();
probj.genMpsf();
probj.phaseretrieve();
probj.genZKresult();
probj.findOTFparam();

probj.genPRfigs('PSF');
set(gcf,'PaperPositionMode','auto')
ImgName='-PSF';
print(gcf, '-r300', '-djpeg', [FileDir,FileName(1:end-4),ImgName])

probj.genPRfigs('pupil');
probj.genPRfigs('zernike');


