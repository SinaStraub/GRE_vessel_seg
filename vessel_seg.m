function [vein_seg]=vessel_seg(magh,phh,qsmh,maskh,echoes,scales,vein_art,shearSys,do_r2,r2,voxelsize,padsize,do_field,mask01,it)

% Compute a vein (or arteries) segmentation (output=vein_seg) from multi-echo or single-echo
% gradient echo data. As the algorithm is implemented in a modular way, it
% is possible to replace the functions for background suppression, calculation of the 
% vesselness function, the scale-wise representation by different methods 
% Implementation based on (please kindly cite the following paper if you use 
% this program, or any code extended from this program:
% [1] Authors, ``A gradient echo data based vein segmentation algorithm and
% its application for the detection of regional cerebral differences in
% venous susceptibility``, Journal 2021,...
%
% Required input data:
% magh: single or multi-echo gradient echo magnitude images, echoes need to be
% fourth dimension.
% phh: single or multi-echo gradient echo phase images, echoes need to be
% fourth dimension.
% qsmh: susceptibility maps from one or multiple echoes separately or
% combined, echoes need to be fourth dimension.
% maskh: brain mask(s) used for susceptibility map calculation, echoes need
% to be fourth dimension.
% echoes: for multi-echo data the echo number of the data that will be used
% for the vein segmentation
% scales: number of scales used for the shearlet transform
% kernel: kernel size for the local thersholding
% vein_art: 0 generates a vein segmentation, a number ~0 an arteries segmentation
% shearSys: used by make_recon, should be pre-computed with make_Shear_sys
% do_r2: parameter that indicates the availibility of R2* data (1 is available)
% r2: R2* map
% voxelsize:in mm [x y z]
% do_field: 1 use regularization function for vesselness as penalty for
% field inhomogeneites, 0 do not use.
% it: vein segmentation for the lowest scale is iterated it-times to improve the
% completeness of segmented large veins.
%
%  Example:
%  First load all required data (including mag), then compute shearSys by [shearSys]=make_Shear_sys(mag,scales);
%  vein_seg=vessel_seg(mag,ph,qsm,mask_erode,2,4,[41,21,5,5],0,shearSys,1,r2,[0.5 0.5 0.6],[6,6,6],1,mask01,0);
%  Uses the second echo of the provided data, 4 scales with thresholds for
%  local thersholdung of 41, 21, 5,5 pixels to compute a vein segmentation.
%
% Dependencies:
% Functions provided: getfield1, BG_supp, make_recon, make_Shear_sys, PFT3Dmod,
% QSM_SWI, vessel_seg_loc, remove_culster
% Third party code:
% Laplacian phase unwrapping used in getfield1:
% https://people.eecs.berkeley.edu/~chunlei.liu/software.html
% 3D Shearlet transform required in make_Shear_sys and make_recon:
% http://shearlab.math.lmu.de/software#shearlab3D (ShearLab3D v1.1)
% Vesselness function which is the source for PFT3Dmod: https://github.com/Haifafh/MFAT
% 
% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 18.012.2021 V1.1.1
%
%Start 1: Handle multi-echo data
if size(magh,4)>1
    mag=magh(:,:,:,echoes);
else mag=magh;
end
if size(phh,4)>1
    ph=phh(:,:,:,echoes);
else ph=phh;
end
if size(maskh,4)>1
    mask=maskh(:,:,:,1);
else mask=maskh;
end
if size(qsmh,4)>1
    qsm=qsmh(:,:,:,echoes);
else qsm=qsmh;
end
%End 1: Handle multi-echo data
%Start 2: Set default threshold values for local segmentations, brain extract
% input data, initialize segmentation variables
t1=2.5;
t2=2.5;
t3=2.5;
t4=2.4;
t5=2.6;
t6=2.7;
t7=2.8;
t8=2.9;
t9=2;
t10=2.3;
ker_art=round([12.5,9.5,4.5]/voxelsize(1));%replace by mean(voxelsize) by for strongly anisotropic resolutions
kernel=round([21.5,10.5,2.5]/voxelsize(1));%replace by mean(voxelsize) by for strongly anisotropic resolutions
mag=mag.*mask;
qsm=qsm.*mask;
if do_r2==1
    r2=r2.*mask;
end
X=size(mag,1);
Y=size(mag,2);
Z=size(mag,3);
vein_seg=zeros(size(mag));
vein_seg2=zeros(size(mag));
veins_all=zeros(X,Y,Z,scales-1);

% End 2: Set default threshold values for local segmentations, brain extract
% input data, initialize segmentation variables

%Start 3: Compute regularization function for field inhomogeneities
if do_field==1
    [field]=getfield1(ph,voxelsize,padsize,mask01);% corresponds to reg in the publication 2.4
else
    field=ones(size(qsm));
end
%End 3: Compute regularization function for field inhomogeneities

if vein_art==0 %segment veins
    
    %Start 4: Compute tSWI, corresponds to step 2.1 in the publication 
    swi = QSM_SWI(mag, qsm);
    %End 4: Compute tSWI
    
    %Start 5: Background suppression with inverse Hamming filter
    swi_BG=BG_supp(swi);
    qsm_BG=BG_supp(qsm);
    if do_r2==1 %for R2* if available
        r2_BG=BG_supp(r2);
    else
        r2_BG=zeros(size(qsm));
    end
    %End 5: Background suppression with inverse Hamming filter, 
    % corresponds to step 2.2 in the publication 
    
    %Start 6: Normalize
    minn=min(qsm_BG(:));
    maxx=max(qsm_BG(:));
    qsm_BGn=(qsm_BG-minn)./(maxx-minn);
    minn=min(swi_BG(:));
    maxx=max(swi_BG(:));
    swi_BGn=(swi_BG-minn)./(maxx-minn);
    if do_r2==1 %for R2* if available
        minn=min(r2_BG(:));
        maxx=max(r2_BG(:));
        r2_BGn=(r2_BG-minn)./(maxx-minn);
    end
    %End 6: Normalize, corresponds to step 2.2 in the publication 
    
    %Start 7: Compute shearlet representation of normalized and background
    %suppressed data.
    %corresponds to step 2.3 in the publication 
    [recon_qsm]=make_recon(qsm_BGn,shearSys);
    [recon_swi]=make_recon(swi_BGn,shearSys);
    if do_r2==1 %for R2* if available
        [recon_r2]=make_recon(r2_BGn,shearSys);
    else
        recon_r2=0;
    end
    % End 7: Compute shearlet representation of normalized and background
    % suppressed data.
    
    % Start 8: Compute vesselness (fixed parameters are used) and threshold locally
    % corresponds to step 2.4 in the publication 
    [vessels_all1,~] = PFT3Dmod(recon_swi,0.02,0.35,0.3,false,scales-1);
    if do_r2==1 %for R2* if available
        [vessels_all2,~] = PFT3Dmod((recon_r2.*recon_qsm),0.02,0.35,0.3,true,scales-1);
        %
        [vein_seg2_1,~] = vessel_seg_loc( round(20.5/voxelsize(1)), vessels_all2(:,:,:,scales-2).*field,1,t4,t4,t4);
        [vein_seg2_2,~] = vessel_seg_loc(  round(15.5/voxelsize(1)), vessels_all2(:,:,:,scales-1).*field,1,t4,t4,t4);
        [vein_seg2_3,~] = vessel_seg_loc(  round(7.5/voxelsize(1)), vessels_all2(:,:,:,scales-1).*field,1,t7,t7,t7);
        vein_seg2(vein_seg2_1==1 | vein_seg2_2==1 | vein_seg2_3==1)=1;
    else
        [vessels_all2,~] = PFT3Dmod(recon_qsm,0.02,0.35,0.3,true,scales-1);
        [vein_seg2_1,~] = vessel_seg_loc( round(20.5/voxelsize(1)), max(cat(4,vessels_all1(:,:,:,scales-2),vessels_all2(:,:,:,scales-2)),[],4).*field,1,t4,t4,t4);
        [vein_seg2_2,~] = vessel_seg_loc( round(15.5/voxelsize(1)), max(cat(4,vessels_all1(:,:,:,scales-1),vessels_all2(:,:,:,scales-1)),[],4).*field,1,t4,t4,t4);
        [vein_seg2_3,~] = vessel_seg_loc( round(7.5/voxelsize(1)), max(cat(4,vessels_all1(:,:,:,scales-1),vessels_all2(:,:,:,scales-1)),[],4).*field,1,t7,t7,t7);
        vein_seg2(vein_seg2_1==1 | vein_seg2_2==1 | vein_seg2_3==1)=1;
    end
    vein_seg2=vein_seg2.*mask;
    % End 8: Compute vesselness and threshold locally
    
    % Start 9: Local thresholding of different scales of tSWI, QSM and if
    % available R2*
    for l=1:scales-1
        % Start 9a: Initialize variables
        help6=zeros(X,Y,Z);
        remove=zeros(X,Y,Z);
        % End 9a: Initialize variables
        
        % Start 9b: Choose scale
        help2=recon_qsm(:,:,:,l);
        if do_r2==1 %for R2* if available
            helpr2=recon_r2(:,:,:,l);
        end
        
        help12=recon_swi(:,:,:,l);
        hh1=help12(:);
        mm1=mask(:);
        help12(mask==0)=mean(hh1(mm1==0))+0.2*std(hh1(mm1==0));% increase values outside the brain
        % End 9b: Choose scale
        
        % Start 9c: Local thresholding scale-wise
        if l==1 %corresponds to step 2.5a in the publication 
            [vein_seg1,~] = vessel_seg_loc( kernel(l), help2,1,t5,t5,t5);
            
            if do_r2==1 %for R2* if available
                [vein_segr2,~] = vessel_seg_loc( kernel(l), helpr2,1,t6,t6,t6);
                vein_seg_r20=vein_segr2;
                vein_seg_r20=remove_culster(vein_seg_r20,1);
                vein_seg_r20=imfill( vein_seg_r20,'holes');%fill holes
            end      
            
            % Start 9d: Iterate local thresholding for scale 1
            for tt=1:it
                
                help2(vein_seg1==1)=mean(help2(mask==1 & vein_seg1==0));
                [vein_seg1h,~] = vessel_seg_loc( kernel(1), help2,1, t6,t6,t6);
                vein_seg1(vein_seg1h==1)=1;
                
                if do_r2==1 %for R2* if available
                    helpr2(vein_segr2==1)=mean(helpr2(mask==1& vein_segr2==0));
                    [vein_segr2h,~] = vessel_seg_loc( kernel(1), helpr2,1, t8,t8,t8);
                    vein_segr2(vein_segr2h==1)=1;
                end
            end
            % End 9d: Iterate local thresholding for scale 1
            
            vein_seg_qsm0=vein_seg1;
            vein_seg_qsm0=remove_culster(vein_seg_qsm0,1);
            vein_seg_qsm0=imfill( vein_seg_qsm0,'holes');
            
            if do_r2==1 %for R2* if available
                vein_seg_r20=vein_segr2;
                vein_seg_r20=remove_culster(vein_seg_r20,1);
            end
        elseif l==scales-1
            [vein_seg1,~] = vessel_seg_loc( kernel(l), help2,1,t10,t10,t10);
            
            [vein_seg3,~] = vessel_seg_loc( kernel(l), help12,0,t10,t10,t10);

        else
            [vein_seg1,~] = vessel_seg_loc( kernel(l), help2,1,t1,t2,t3);
            
            [vein_seg3,~] = vessel_seg_loc( kernel(l), help12,0,t1,t2,t3);
        end
        % End 9c: Local thresholding scale-wise
        % End 9: Local thresholding of different scales of tSWI, QSM and if
        % available R2*
        
        % Start 10: Segment and calculated segmentation for calcified structures
        % (dark in SWI, dark on QSM)
        %corresponds to step 2.5b in the publication, remove_... corresponds to seg_remove 
        [remove_qsm,~] = vessel_seg_loc( round(30.5/voxelsize(1)), help2,0,t9,t9,t9);
        [remove_swi,~] = vessel_seg_loc( round(30.5/voxelsize(1)), help12,0,t9,t9,t9);
        remove((remove_qsm==1 & remove_swi==1))=1;
        
        remove((remove_qsm==1 & remove_swi==1))=1;
        % End 10: Segment and calculated segmentation for calcified structures
        % (dark in SWI, dark on QSM
        
        % Start 11: Include/ remove voxels in/ from segmentation
        if l==1
            help6((vein_seg1==1 & vein_seg2==1 & remove==0) )=1;
        else
            help6((vein_seg1==1 & vein_seg2==1 & remove==0) | (vein_seg2==1 & vein_seg3==1 & remove==0) )=1;
        end
        
        veins_all(:,:,:,l)=help6;
        vein_seg(help6==1)=1;
        % Start 11a: Include voxels for lowest scale
        if do_r2==1 %for R2* if available
            vein_seg( vein_seg_qsm0==1 & vein_seg_r20==1)=1;
        else
            vein_seg(vein_seg_qsm0==1)=1;
        end
        % End 11a: Include voxels for lowest scale
    end
    vein_seg(mask==0)=0;
    vein_seg=remove_culster(vein_seg,1);% remove clusters smaller than 1 voxel
    % End 11: Include/ remove voxels in/ from segmentation
    
else% segment arteries
    % Start: Normalize
    minn=min(mag(:));
    maxx=max(mag(:));
    magn=(mag-minn)./(maxx-minn);
    % End: Normalize
    
    % Start: Compute shearlet representation
    [recon_swi]=make_recon(magn,shearSys);
    % End: Compute shearlet representation
    
    %Start: Vesselness and local thresholding
    [vessels_all1,~] = PFT3Dmod(recon_swi,0.02,0.35,0.3,true,scales-1);
    [ ~,vein_seg2_1] = vessel_seg_loc( ker_art(1), vessels_all1(:,:,:,1).*(1-field),0,t7,t8,t9);
    [ ~,vein_seg2_2] = vessel_seg_loc( ker_art(2), vessels_all1(:,:,:,2).*(1-field),0,t7,t8,t9);
    [ ~,vein_seg2_3] = vessel_seg_loc( ker_art(3), vessels_all1(:,:,:,3).*(1-field),0,t7,t8,t9);
    vein_seg2(vein_seg2_1==1 | vein_seg2_2==1 |  vein_seg2_3==1)=1;
    
    vein_seg2=vein_seg2.*mask;
    %End: Vesselness and local thresholding
    
    % Start: Local threseholding of different scales
    for l=1:scales-1
        help12=recon_swi(:,:,:,l);
        hh1=help12(:);
        mm1=mask(:);
        help12(mask==0)=mean(hh1(mm1==0))+2*std(hh1(mm1==0));
        
        [~,vein_seg3] = vessel_seg_loc( kernel(l), help12,0,t7,t8,t9);
        vein_seg( vein_seg2==1 & vein_seg3==1  )=1;
    end
    % End: Local threseholding of different scales
    
    % Start: Remove clusters smaller than two voxels, and brain extract
    vein_seg=remove_culster(vein_seg,2);
    vein_seg=remove_culster(vein_seg,1);
    vein_seg(mask==0)=0;
    % End: Remove clusters smaller than two voxels, and brain extract
end




