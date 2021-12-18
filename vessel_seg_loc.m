function [seg_vein,seg_artery] = vessel_seg_loc( kernel_size, image,contrast,t1,t2,t3)
%contrast==0 -> mag, swi etc. where veins are dark, arteies bright
%contrast==1 -> qsm etc. where veins are bright, arteies darker

% Examples:[seg_vein_loc,seg_artery_loc] = vessel_seg_loc( 11, mag,0,5,5,5),
%[seg_vein_loc,~] = vessel_seg_loc( 11, qsm,0,5,5,5)

% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de

%kernelsize needs to be odd, so that kernel is centered at current element
seg_vein=zeros(size(image));
seg_artery=zeros(size(image));

if nargin<5
    t1=2.5;
    t2=2.5;
    t3=2.5;
end

smatrix=size(image);

[Y,X,Z]=meshgrid(-smatrix(2)/2:smatrix(2)/2-1,...
                 -smatrix(1)/2:smatrix(1)/2-1,...
                 -smatrix(3)/2:smatrix(3)/2-1);
%kernel mean
kernel=zeros(smatrix);
kernel (( abs(X).^2  +abs(Y).^2 +abs(Z).^2 )<=kernel_size^2)=1; 
kernel_help=sum(kernel(:));
kernel=kernel./kernel_help;
fkernel=fftn(fftshift(kernel));
fimage=fftn(fftshift(image));
filtered=fkernel.*fimage;
help_mean=ifftshift(ifftn(filtered));
stdh1=(image-help_mean).^2;
fstd=fftn(fftshift(stdh1));
stdh2=fkernel.*fstd;
help_std=sqrt(ifftshift(ifftn(stdh2)));

if contrast==0

    seg_vein(image<help_mean-t1*help_std)=1; 
    seg_artery(image>help_mean+t2*help_std)=1; 

elseif contrast==1
  seg_vein(image>help_mean+t3*help_std)=1;   
end

end
