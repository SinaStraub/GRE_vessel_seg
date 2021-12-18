% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 07.10.2021 V1.1

addpath(genpath('path_for_vessel_seg_code_and_thid_party_code'))% adjust to fit your path! All dependencies must be included.
%Example 1: 1.0 x 1.0 x 1.2 mm³ resolution, TE=6.28
load('data_ex1') 
[shearSys]=make_Shear_sys(mag,4);
[vein_seg]=vessel_seg(mag,ph,qsm,mask,2,4, 0,shearSys,1,r2,[1 1 1.2],[6 6 6],1, mask01,20);

%save result before proceeding

%Example 2: 0.5 x 0.5 x 1.8 mm³ resolution, combined TE from 12 echoes: TE=3.14/ 6.28/ 9.42/ 11.99/ 14.56/ 17.13/ 19.7/ 22.27/ 24.84/ 27.41/ 29.98/ 32.55
load('data_ex2') 
[shearSys]=make_Shear_sys(mag(:,:,:,1),4);
[vein_seg]=vessel_seg(mag,ph,qsm,mask,1,4,0,shearSys,1,r2,[0.5 0.5 1.8],[6 6 6],1, mask01,20);

   
