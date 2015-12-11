function [tfced] = tfce_transform(img,H,E,C)
%TFCE_TRANSFORM performs threshold free cluster enhancement on image
%   [tfced] = tfce_transform(img,H,E,C) performs a threshold free cluster
%   enhancement on the image 'img' as per Smith & Nichols (2009).
%   -- H height exponent, default = 2
%   -- E extent exponent, default = 0.5
%   -- C connectivity, default = 18 (6 = faces, 18 = sides, 26 = corners)
%   The H & E default parameter settings match FSL's randomise/fslmaths.
%   The C default setting matches SPM's default cluster definition. To
%   match FSL's default, use 6 instead.


end

