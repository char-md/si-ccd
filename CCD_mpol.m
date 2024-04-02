function imOut = CCD_mpol(imIn1,imIn2,winSize)

% This function implements the multi-polarization based change detection
% method described in 
% Change Detection for Multi-Polarization Multi-Pass SAR
% Algorithms for Synthetic Aperture Radar Imagery XII
% Novak, Leslie M.
% vol. 5808, pp. 234-246, SPIE, May 2005
% doi: 10.1117/12.609894
%


%% Initialization
% Image size
% Note the inputs are image tensors
nx = size(imIn1,1);
ny = size(imIn1,2);

% No. of channels or polarizations in the inputs
ncomp = size(imIn1,3);

% Initialize output image
imOut = ones(nx,ny);

% Temporary variables
f = zeros(ncomp, (2*winSize+1)^2);
g = zeros(ncomp, (2*winSize+1)^2);


% Threshold
thr = 1e-3;


%% Multi-polarization change map
% TODO: Figure out how to make this faster
% Proceed pixel-by-pixel
for ix = winSize+1:nx-winSize
    for iy = winSize+1:ny-winSize
    
        % Collapse the channels into a single vector for both input images
        for icomp = 1:ncomp
            f(icomp,:) = reshape(imIn1(ix-winSize:ix+winSize,iy-winSize:iy+winSize,icomp), 1, (2*winSize+1)^2 );
            g(icomp,:) = reshape(imIn2(ix-winSize:ix+winSize,iy-winSize:iy+winSize,icomp), 1, (2*winSize+1)^2 );
        end

        % Compute outer products
        % This implements a variant of (13) in the above manuscript
        time1 = f*f'/((2*winSize+1)^2);        
        time2 = g*g'/((2*winSize+1)^2);

        % Check for small denominators and mask
        den = ((det(0.5*(time1+time2)))^2);
        if (abs(den)>thr)
            imOut(ix,iy) = ( det(time1)*det(time2) )/den;
        end
    end
end

% There can be complex entries in the result; return only the magnitude
imOut = abs(imOut);

return