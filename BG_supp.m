function [new_image] = BG_supp(mag)
% function computes inverse Hamming filtered image (new_image) from the
% input 3D image mag, the implementation is based on the publication:
% Jin, Z.Y., Xia, L., Zhang, M.M., Du, Y.P.P., 2014. Background-Suppressed 
% MR Venography of the Brain Using Magnitude Data: A High-Pass Filtering Approach.
% Computational and Mathematical Methods in Medicine.

% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 27.03.2021 V1.1

image=mag+1.i*0;
fimage=ifftshift(fftn(fftshift(image)));


gridx=repmat([-size(mag,1)/2:-1,1:size(mag,1)/2]',1,size(mag,2),size(mag,3));
gridy=repmat([-size(mag,2)/2:-1,1:size(mag,2)/2],size(mag,1),1,size(mag,3));
help(1,1,:)=[-size(mag,3)/2:-1,1:size(mag,3)/2];
gridz=repmat(help,size(mag,1),size(mag,2),1);
Hx=80;
Hy=80;
Hz=80;


filter=ones(size(mag));
help2=(gridx.^2/Hx^2+gridy.^2/Hy^2+gridz.^2/Hz^2);
filter((gridx.^2/Hx^2+gridy.^2/Hy^2+gridz.^2/Hz^2)<=1)= 0.6*(1-cos(pi*sqrt(help2((gridx.^2/Hx^2+gridy.^2/Hy^2+gridz.^2/Hz^2)<=1))));
filtered=fimage.*filter;
new_image=real(fftshift(ifftn(ifftshift(filtered))));

end