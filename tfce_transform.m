function [tfced] = tfce_transform(img,H,E,C,ndh)
%TFCE_TRANSFORM performs threshold free cluster enhancement on image
%   [tfced] = tfce_transform(img,H,E,C) performs a threshold free cluster
%   enhancement on the image 'img' as per Smith & Nichols (2009).
%   -- H height exponent, default = 2
%   -- E extent exponent, default = 0.5
%   -- C connectivity, default = 6 (6 = surface, 18 = edge, 26 = corner)
%   -- ndh step number for cluster formation, default = 100
%   More steps will be more precise but will require more time and memory.
%   The H & E default parameter settings match FSL's randomise/fslmaths.
%   The C default setting matches FSL's ramdomise default setting. To
%   match SPM's default cluster forming, use 18 instead.

% set cluster thresholds
threshs = linspace(0,max(img(:)),ndh+2);
threshs = threshs(2:(ndh+1));
dh = diff(threshs);
dh = dh(1);

% find positive voxels (greater than first threshold)
nvox = length(img(:));
posvoxels = (img>=threshs(1));
posinds = find(posvoxels);
pnum = length(posinds);

% find connected components
clusts = NaN(pnum,ndh);
for h = 1:ndh
    cc = bwconncomp(img>=threshs(h),C);
    curclust = NaN(nvox,1);
    for c = 1:cc.NumObjects
        curclust(cc.PixelIdxList{c}) = c;
        clusts(:,h) = curclust(posinds);
    end
end

% run through positive voxels
tfced = img;
posvals = img(posinds);
for p = 1:pnum
    nsupthr = sum(posvals(p)>threshs);
    thrvals = NaN(nsupthr,1);
    for nh = 1:nsupthr
        cid = clusts(p,nh);
        extent = sum(clusts(:,nh)==cid);
        thrvals(nh) = (extent^E)*(dh^H);
    end
    tfced(posinds(p)) = sum(thrvals);
end

end

