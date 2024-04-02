close all; clear all; clc

%%%%%%%%%
% TODO: Lots of changes to be done
%%%%%%%%%

% Parts of this code are based on the Matlab script files included as part 
% of the GOTCHA CCD challenge problem dataset
%
% The original version of the scripts were written by 
% Steven Scarborough and LeRoy Gorham (AFRL/RYAP)
% Email:  steven.scarborough@wpafb.af.mil / leroy.gorham@wpafb.af.mil


% Define the path to the base directory of the dataset
datadir = '~/Downloads/ChangeDetectionDataset_GOTCHA/CCD-CP-XPol/HH/';

% Define input pass names here
pass1 = 'FP0120';             % Which pass is the reference image
pass2 = 'FP0124';             % Which pass is the mission image

% Define change detection parameters here
CCDwinSize = 3;               % Window size for CCD (see help file for CCDbasic.m)

% Define image parameters here
minRow = 1300;                % Minimum row value
maxRow = 3300;                % Maximum row value
minCol = 1200;                % Minimum column value
maxCol = 3200;                % Maximum column value


dyn_range = 60;               % dB of dynamic range to display on SAR image outputs



% Determine the file names of the input files
% Reads in HV polarization
im1filename = sprintf('%s%s/c00007a283p50.mat',datadir,pass1);
im2filename = sprintf('%s%s/c00007a283p50.mat',datadir,pass2);


% Load in the reference image
data1 = load(im1filename);

% Crop the image to requested size
im1 = data1.SARdataOut(minRow:maxRow,minCol:maxCol);

% Display the reference image
figure
imagesc([],[],20*log10(abs(im1)./max(abs(im1(:)))),[-dyn_range 0]);
axis image
title('Reference Image');
colorbar


%% Edges in the reference image
% First some filtering to remove speckle
% This is an edge-preserving anisotropic diffusion based filtering
% implementation
im1_filt = specklefilt(abs(im1),DegreeOfSmoothing=0.6,NumIterations=50);

% Edge detection
im1_edge = sqrt(ComputeEdgesImage(im1_filt,'xdim',CCDwinSize).^2 + ComputeEdgesImage(im1_filt,'ydim',CCDwinSize).^2 );


% Display edges in the reference image
figure
imagesc(im1_edge)
axis image
title('Reference Image - Edges');
colorbar


%% 

% Load in the mission image
data2 = load(im2filename);

% Crop the image to requested size
im2 = data2.SARdataOut(minRow:maxRow,minCol:maxCol);

% Display the mission image
figure
imagesc([],[],20*log10(abs(im2)./max(abs(im2(:)))),[-dyn_range 0]);
axis image
title('Mission Image');
colorbar


%% Edges in mission image
% First some filtering to remove speckle
% This is an edge-preserving anisotropic diffusion based filtering
% implementation
im2_filt = specklefilt(abs(im2),DegreeOfSmoothing=0.6,NumIterations=50);

% Edge detection
im2_edge = sqrt( ComputeEdgesImage(im2_filt,'xdim',CCDwinSize).^2 + ComputeEdgesImage(im2_filt,'ydim',CCDwinSize).^2 );


% Display edges in the reference image
figure
imagesc(im2_edge)
axis image
title('Mission Image - Edges');
colorbar


drawnow

%%


% Perform CCD using CCDbasic.m
CCDimage = CCDbasic(im1,im2,CCDwinSize,'MLE',true);




% Enhanced CCD

% Low-RCS masking
% Apply entropy filtering
mask_fnc_im1 = entropyfilt(abs(im1),true(2*CCDwinSize+1));
mask_fnc_im2 = entropyfilt(abs(im2),true(2*CCDwinSize+1));
% the second argument is the neighborhood size to use

% Reconcile results from the two images
mask_fnc = min(mask_fnc_im1, mask_fnc_im2);

% Rescale to [0,1]
mask_fnc = 1-rescale(mask_fnc);

    
% Masking threshold
% TODO: make this an optional parameter input with a default value
th = 0.15;    

% Find regions which satisfy the masking condition
mask_flag = (mask_fnc<=th);

% Enforce such regions to have unit signal value (so that they don't
% register in the change map)
im1(mask_flag) = 1;
im2(mask_flag) = 1;



% Multipolar CCD with edge information
CCDimage_mpol = CCD_mpol(cat(3,im1,im1_edge, sign(im1)), cat(3,im2,im2_edge, sign(im2)), CCDwinSize);


% Image connectivity analysis for false positive mitigation
% Initialize connectivity matrix
ncon = zeros(size(im1));

% Offset for speedup
offset = 5;

% Break up change map into regions (square patches)
reg_size = 50;         % size of patch

for ix = reg_size/2+1:offset:size(im1,1)-reg_size/2
    for iy = reg_size/2+1:offset:size(im1,2)-reg_size/2
        % Extract image
        im = CCDimage_mpol(ix-reg_size/2:ix+reg_size/2, iy-reg_size/2:iy+reg_size/2);
        
        % Connectivity
        % tmp = bwconncomp(im2bw(im));
        tmp = bwconncomp(imbinarize(im, 'adaptive', 'ForegroundPolarity', 'bright'));
        ncon(ix-reg_size/2:ix+reg_size/2, iy-reg_size/2:iy+reg_size/2) = tmp.NumObjects;
    end
end

% Connectivity threshold
th = 0.35*max(ncon(:));

% Generate mask using connectivity information
mask = (ncon<=th);

% Generate post-processed change map
CCDimage_mpol_pp = CCDimage_mpol;
CCDimage_mpol_pp(mask==0) = 1;



% Display the CCD results
figure
imagesc(CCDimage,[0 1]);
axis image
title('Coherent Change Detection');
colormap gray
colorbar


figure
imagesc(CCDimage_mpol_pp,[0 1]);
axis image
title('Coherent Change Detection (w edge information; false positive mitigation)');
colormap gray
colorbar




% Combine the MLE-based and multipolar change maps

% Image connectivity analysis for false positive mitigation of the
% MLE-based change map
% Initialize connectivity matrix
ncon = zeros(size(im1));

% Offset for speedup
offset = 5;

% Break up change map into regions (square patches)
reg_size = 50;         % size of patch

for ix = reg_size/2+1:offset:size(im1,1)-reg_size/2
    for iy = reg_size/2+1:offset:size(im1,2)-reg_size/2
        % Extract image
        im = CCDimage(ix-reg_size/2:ix+reg_size/2, iy-reg_size/2:iy+reg_size/2);
        
        % Connectivity
        % tmp = bwconncomp(im2bw(im));
        tmp = bwconncomp(imbinarize(im, 'adaptive', 'ForegroundPolarity', 'bright'));
        ncon(ix-reg_size/2:ix+reg_size/2, iy-reg_size/2:iy+reg_size/2) = tmp.NumObjects;
    end
end

% Connectivity threshold
th = 0.35*max(ncon(:));

% Generate mask using connectivity information
mask = (ncon<=th);

% Generate post-processed change map
CCDimage_pp = CCDimage;
CCDimage_pp(mask==0) = 1;


% Take a convex combination of the two
CCDimage_comb = 0.75*CCDimage_pp + 0.25*CCDimage_mpol_pp;

% Plot this combined change map image
figure
imagesc(CCDimage_comb,[0 1]);
axis image
title('Coherent Change Detection (composite)');
colormap gray
colorbar
