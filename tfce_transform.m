function [tfced] = tfce_transform(img,varargin)
%TFCE_TRANSFORM performs threshold free cluster enhancement on image
%   [tfced] = tfce_transform(img,H,E,C,ndh) performs threshold free cluster
%   enhancement on the image 'img' as per Smith & Nichols (2009).
%   -- img the 3D image to be transformed
%   -- H height exponent, default = 2
%   -- E extent exponent, default = 0.5
%   -- C connectivity, default = 6 (6 = surface, 18 = edge, 26 = corner)
%   -- ndh step number for cluster formation, default = 100
%   More steps will be more precise but will require more time and memory.
%   The H & E default parameter settings match FSL's randomise/fslmaths.
%   The C default setting matches FSL's ramdomise default setting. To
%   match SPM's default cluster forming, use 18 instead. The transformed
%   image is returned as 'tfced'. Note that for most purposes, the wrapper
%   function tfce_transform_twotailed should be applied to the image to
%   ensure that results are not biased in the positive direction.

% setting defaults
H = 2;
E = .5;
C = 6;
ndh = 100;

% overriding defaults as necessary
if nargin > 1
    H = varargin{1};
end
if nargin > 2
    E = varargin{2};
end
if nargin > 3
    C = varargin{3};
end
if nargin > 4
    ndh = varargin{4};
end

% set cluster thresholds
threshs = linspace(0,max(img(:)),ndh+2);
threshs = threshs(2:(ndh+1));
dh = diff(threshs);
dh = dh(1);

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
vals = sum((clustsize.^E).*(dh^H),2);
tfced = NaN(size(img));
tfced(:) = vals;
tfced(img<0) = img(img<0);

end

