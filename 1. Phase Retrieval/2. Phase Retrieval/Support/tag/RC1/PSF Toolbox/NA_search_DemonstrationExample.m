% NA_search_testDemonstrationExample

% load the Phase Retrieval data (analyzed by PR code)
clear all
a=NA_search_test() % define the object
% set the properties
a.Bpp=0.096; % micron 
a.FOV=128;   % # of pixels (FOV is defined based on the size of the image)
a.Lambda_emm=690;
load('Y:\Farzin\Sequential Imaging Data\2015_7_6_CRLBbeads\PR_result\darkredbead-2015-7-6-16-5-46_PR_optim_single.mat')
a.RawData=obj.Data; % left: defined object above.name of the property in class / right: obj from PR results.Data
a.findOTF_Radius(a.Bpp,a.FOV) % finds radius of OTF of the data and then gives back NA_data
a.showResults(a.FOV,a.Bpp,a.Lambda_emm)

% Note: if you have functions under @folder, don't put them under methods.
% Note: if a function is big, don't put it under methods, but under @folder

