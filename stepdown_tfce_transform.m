function [tfced] = stepdown_tfce_transform(img,H,E,C,ndh)
%STEPDOWN_TFCE_TRANSFORM performs threshold free cluster enhancement
%   [tfced] = stepdown_tfce_transform(img,H,E,C,ndh) performs threshold
%   free cluster enhancement on 'img' as per Smith & Nichols (2009).
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- ndh step number for cluster formation

% set cluster thresholds
threshs = linspace(0,max(img(:)),ndh+2);
threshs = threshs(2:(ndh+1));

% find positive voxels (greater than first threshold)
nvox = length(img(:));

% find connected components
clustsize = zeros(nvox,ndh);
for h = 1:ndh
    cc = bwconncomp(img>=threshs(h),C);
    voxpercc = cellfun(@numel,cc.PixelIdxList);
    for c = 1:cc.NumObjects
        clustsize(cc.PixelIdxList{c},h) = voxpercc(c);
    end
end

% calculate TFCE values and insert
vals = sum((clustsize.^E).*(repmat(threshs,nvox,1).^H),2);
tfced = NaN(size(img));
tfced(:) = vals;

end

