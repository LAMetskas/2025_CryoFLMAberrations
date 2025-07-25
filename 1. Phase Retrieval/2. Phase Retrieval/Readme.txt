

run program 'PR_est_with_true_image_iters' --> phase retrieval for multiple iterations

Input: ims_Ztrue_plane1_s and ims_Ztrue_plane1_n (you already generate in file step 2)

output: measured and phase retrieved PSF
        phase retrieved and Zernike fitted pupil function
        phase_Zernike coefficinet 
        
First, you need change Zpos_plane1 according the real data. In this case Z range is about [-1,1] um.
You can aslo change to [-1.5, 1.5] for new data.


You can adjust all the following parameters based on your own case:
CCD camera properties; Optical system configuration; Phase retrieval and PSF settings; Gaussian blur parameter; max_iterations;

