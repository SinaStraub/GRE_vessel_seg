function [tswi] = QSM_SWI(mag, qsm)
%calculates tSWI from qsm and GRE magnitude 

% This implementation is based:
% Liu, S.F., Mok, K., Neelavalli, J., Cheng, Y.C.N., Tang, J., Ye, Y.Q., Haacke,
% E.M., 2014. Improved MR Venography Using Quantitative Susceptibility-Weighted
% Imaging. Journal of Magnetic Resonance Imaging 40, 698-708.

% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 27.03.2021 V1.1
n=4; 
%tissue thresholds:
thres1= 0.0005;
thres2= 0.45;
%phase/qsm mask:
mask=zeros(size(qsm,1),size(qsm,2),size(qsm,3));
mask(qsm<thres1)=1;
help_mask=ones(size(qsm,1),size(qsm,2),size(qsm,3))-((qsm-thres1)/(thres2-thres1));
mask(qsm>=thres1 & qsm<=thres2)=help_mask(qsm>=thres1 & qsm<=thres2);
tswi=mag.*mask.^n;
end

