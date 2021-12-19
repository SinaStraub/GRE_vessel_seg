GRE_vessel_seg computes a vein (or arteries) segmentation (output=vein_seg) from multi-echo or single-echo
gradient echo data in Matlab. As the algorithm is implemented in a modular way, it
is possible to replace the functions for background suppression, calculation of the 
vesselness function, the scale-wise representation by different methods. 
Implementation based on (please kindly cite the following paper if you use 
this program, or any code extended from this program):

[1] Authors, ``A gradient echo data based vein segmentation algorithm and
its application for the detection of regional cerebral differences in
venous susceptibility``, Journal 2021,...

Required input data:

magh: single or multi-echo gradient echo magnitude images, echoes need to be fourth dimension.

phh: single or multi-echo gradient echo phase images, echoes need to be fourth dimension.

qsmh: susceptibility maps from one or multiple echoes separately or combined, echoes need to be fourth dimension.

maskh: brain mask(s) used for susceptibility map calculation, echoes need to be fourth dimension.

echoes: for multi-echo data the echo number of the data that will be used for the vein segmentation

scales: number of scales used for the shearlet transform

kernel: kernel size for the local thersholding

vein_art: 0 generates a vein segmentation, a number ~0 an arteries segmentation

shearSys: used by make_recon, should be pre-computed with make_Shear_sys

do_r2: parameter that indicates the availibility of R2* data (1 is available)

r2: R2* map

voxelsize:in mm [x y z]

do_field: 1 use regularization function for vesselness as penalty for field inhomogeneites, 0 do not use.

it: vein segmentation for the lowest scale is iterated it-times to improve thecompleteness of segmented large veins.

Example:
First load all required data (including mag) (example data can be downloaded here: https://doi.org/10.5281/zenodo.5791233), then compute shearSys by [shearSys]=make_Shear_sys(mag,scales);

vein_seg=vessel_seg(mag,ph,qsm,mask_erode,2,4,[41,21,5,5],0,shearSys,1,r2,[0.5 0.5 0.6],[6,6,6],1,mask01,0);

Uses the second echo of the provided data, 4 scales with thresholds for local thersholdung of 41, 21, 5,5 pixels to compute a vein segmentation.

Dependencies:
Functions provided: getfield1, BG_supp, make_recon, make_Shear_sys, PFT3Dmod, QSM_SWI, vessel_seg_loc, remove_culster

Third party code:
Laplacian phase unwrapping used in getfield1: https://people.eecs.berkeley.edu/~chunlei.liu/software.html

3D Shearlet transform required in make_Shear_sys and make_recon: http://shearlab.math.lmu.de/software#shearlab3D (ShearLab3D v1.1)

Vesselness function which is the source for PFT3Dmod: https://github.com/Haifafh/MFAT
 
Author: Sina Straub

Email: sina.straub@gmail.com, sina.straub@dkfz.de

Date: 18.12.2021 V1.1.1
