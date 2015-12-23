function [tfced] = stepdown_tfce_transform_twotailed(img,H,E,C,ndh)
%STEPDOWN_TFCE_TRANSFORM_TWOTAILED convenience wrapper for tfce_transform
%   Performs threshold free cluster enhancement of both positive and
%   negative sides of an image.
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- ndh step number for cluster formation

% call tfce on both sides of the image
tfced_pos = stepdown_tfce_transform(img,H,E,C,ndh);
tfced_neg = stepdown_tfce_transform(-img,H,E,C,ndh);
tfced = tfced_pos-tfced_neg;

end

