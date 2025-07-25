classdef NA_search_test < handle
    % This class is being used to find Na/lambda (= radius of pupil function)
    % Help: This code takes in the Phase Retrieval results and gives back the 
    
    properties
        RawData;       % Raw bead data 
        ShowSteps = 1;
        Bpp;  %Back-projected pixel size (micron)
        OTF_Radius;
        FOV;
        Lambda_emm;
        Limit1;  % the fitting limit for OTF of the emitter at the infocus plane, in order to find NA
        Limit2;
    end
    
    %Note: The NA found in experiment is 1.4946
    
    methods  % write a constructor also
                    
    end
    
end

