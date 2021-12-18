function  [recon]=make_recon(swi,shearletSystem)
% This implementation is based on code from
% ShearLad3Dv11\Examples\3D\SLExampleVideoDenoising (http://shearlab.math.lmu.de/software) based on the
% publication: G. Kutyniok, W.-Q. Lim, R. Reisenhofer
% ShearLab 3D: Faithful Digital SHearlet Transforms Based on Compactly Supported Shearlets.
% ACM Trans. Math. Software 42 (2016), Article No.: 5.

% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 27.03.2021 V1.1
sizex=size(swi,1);
sizey=size(swi,2);
sizez=size(swi,3);

levels=[1,1,1,1];

l=length(levels);

coeffs = SLsheardec3D(swi,shearletSystem);

for k=1:l-1
help_coeff=zeros(size(coeffs,1),size(coeffs,2),size(coeffs,3),length(find(shearletSystem.shearletIdxs(:,2)==k)));
help_coeff(:,:,:,shearletSystem.shearletIdxs(:,2)==k)=coeffs(:,:,:,shearletSystem.shearletIdxs(:,2)==k);
thresholdingFactor(k)=mean(help_coeff(:))+0.2*std(help_coeff(:));
clear help_coeff
end
sigma=[30,20,10];

recon=zeros(sizex,sizey,sizez,l-1);
coeff_n=coeffs;

for k=1:l-1
    help=shearletSystem.shearletIdxs(:,2)>k+1;
    help1=sum(help);
    
    coeff_n(:,:,:,shearletSystem.shearletIdxs(:,2)>k+1)=zeros(sizex,sizey,sizez,help1);
    %%thresholding
    coeff_n = coeff_n.*(abs(coeff_n) > thresholdingFactor(k)*reshape(repmat(shearletSystem.RMS,[sizex*sizey*sizez 1]),[sizex,sizey,sizez,length(shearletSystem.RMS)])*sigma(k));

    recon(:,:,:,k)=SLshearrec3D(coeff_n, shearletSystem);
    coeff_n=coeffs;
end

end

