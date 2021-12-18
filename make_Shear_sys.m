function  [shearletSystem]=make_Shear_sys(swih,nScales)
% This implementation is based on code from
% ShearLad3Dv11\Examples\3D\SLExampleVideoDenoising (http://shearlab.math.lmu.de/software) based on the
% publication: G. Kutyniok, W.-Q. Lim, R. Reisenhofer
% ShearLab 3D: Faithful Digital SHearlet Transforms Based on Compactly Supported Shearlets.
% ACM Trans. Math. Software 42 (2016), Article No.: 5.

% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 27.03.2021 V1.1
if size(swih,4)>1
    swi=swih(:,:,:,1);
else swi=swih;
end

sizex=size(swi,1);
sizey=size(swi,2);
sizez=size(swi,3);

levels=ones(1, nScales);
useGPU=0;

directionalFilter = modulate2(dfilters('cd','d')./sqrt(2),'c');
shearletSystem = SLgetShearletSystem3D(useGPU, sizex, sizey, sizez, nScales, levels,0,directionalFilter);

end