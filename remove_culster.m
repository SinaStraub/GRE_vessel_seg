function vein_seg=remove_culster(vein_seg,c) 
% removes clusters smaller than c voxels using the Matlab function
% bwconncomp

% Author: Sina Straub
% Email: sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 27.03.2021 V1.1

CC = bwconncomp(vein_seg);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    TT=sum(numPixels(numPixels ==c) );
    for t=1:TT
    CC = bwconncomp(vein_seg);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = min(numPixels);
    vein_seg(CC.PixelIdxList{idx}) = 0;
    end
    
end