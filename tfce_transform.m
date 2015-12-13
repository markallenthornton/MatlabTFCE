function [tfced] = tfce_transform(img,H,E,C,ndh)
%TFCE_TRANSFORM performs threshold free cluster enhancement on image
%   [tfced] = tfce_transform(img,H,E,C) performs a threshold free cluster
%   enhancement on the image 'img' as per Smith & Nichols (2009).
%   -- H height exponent, default = 2
%   -- E extent exponent, default = 0.5
%   -- C connectivity, default = 6 (6 = surface, 18 = edge, 26 = corner)
%   -- ndh step number for cluster formation; more steps will be more
%   precise but slower to compute.
%   The H & E default parameter settings match FSL's randomise/fslmaths.
%   The C default setting matches FSL's ramdomise default setting. To
%   match SPM's default cluster forming, use 18 instead.

%% identify occupied voxels and find indices
posvoxels = (img>0);
posinds = find(posvoxels);
bsize = size(brainvoxels);
pnum = length(posinds);
nvox = length(img(:));

%% generate neighborhood 
neighborhood = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
if C > 6
    neighborhood = [neighborhood;
                    1 1 0; 1 -1 0; 1 0 1; 1 0 -1;
                    -1 1 0; -1 -1 0; -1 0 1; -1 0 -1;
                    0 1 1; 0 1 -1; 0 -1 1; 0 -1 -1]; 
end
if C > 18
        neighborhood = [neighborhood;
                    1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1;
                    -1 -1 -1; -1 1 -1; -1 -1 1; -1 1 1]; 
end

%% finding clusters
% set cluster thresholds
threshs = linspace(0,max(img),ndh+1);
threshs = threshs(1:ndh);

pxl = {};
for h = 1:ndh
    cc = bwconncomp(img,C);
    for c = 1:cc.NumObjects
        
    end
end

% use bwconncomp to segment the image multiple times at different positive
% thresholds (e.g. 0:.1:max(img))

%% run through voxels
for v = 1:vnum
    
    % define voxel of interest and extract value
    vcent = brainxyz(v,:);
    vval = img(vcent);
    
    % define immediate neighbor indices
    neigbhors = repmat(brainxyz,C,1) + neighborhood;
    ninds = ind2sub(bsize,neighbors);

    % obtain discrete neighbor thresholds
    nthreshs = img(ninds);
    nthreshs = [nthreshs(nthreshs>0 & nthreshs<vval) vval];
    
    % step through adjacent, building up multicluster on the weigh
end

end

