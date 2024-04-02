close all; clear; clc

% This code is based on the Matlab script files included as part of the 
% GOTCHA CCD challenge problem dataset
%
% The original version of the scripts were written by 
% Steven Scarborough and LeRoy Gorham (AFRL/RYAP)           %
% Email:  steven.scarborough@wpafb.af.mil / leroy.gorham@wpafb.af.mil

% for pretty pictures
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex');


%% Problem Parameters
% Define the path to the base directory of the dataset
% datadir = '../HH/';
datadir = '~/Downloads/ChangeDetectionDataset_GOTCHA/CCD-CP-XPol/HH/';

% Define input pass names here
pass1 = 'FP0120';             % Which pass is the reference image
pass2 = 'FP0124';             % Which pass is the mission image

% Define change detection parameters here
CCDwinSize = 03;              % Window size for CCD (see help file for CCDbasic.m)

% Define image parameters here
minRow = 1300;                % Minimum row value
maxRow = 3300;                % Maximum row value
minCol = 1200;                % Minimum column value
maxCol = 3200;                % Maximum column value

dyn_range = 60;               % dB of dynamic range to display on SAR image outputs


%% Read in and display data 
% Determine the file names of the input files
im1filename = sprintf('%s%s/c00007a283p50.mat',datadir,pass1);
im2filename = sprintf('%s%s/c00007a283p50.mat',datadir,pass2);

% Load in the reference image
data1 = load(im1filename);

% Crop the image to requested size
im1 = data1.SARdataOut(minRow:maxRow,minCol:maxCol);

% Display the reference image
figure(1)
imagesc([],[],20*log10(abs(im1)./max(abs(im1(:)))),[-dyn_range 0]);
axis image
title('Reference Image');
colorbar

% Load in the mission image
data2 = load(im2filename);

% Crop the image to requested size
im2 = data2.SARdataOut(minRow:maxRow,minCol:maxCol);

% Display the mission image
figure(2)
imagesc([],[],20*log10(abs(im2)./max(abs(im2(:)))),[-dyn_range 0]);
axis image
title('Mission Image');
colorbar


%% Compute change map without low-RCS masking
% Perform CCD using the MLE method
CCDimage = CCDbasic(im1,im2,CCDwinSize,'MLE',false);
% the last argument determines whether to implement low-RCS masking

% Display the CCD results
figure;
imagesc(CCDimage,[0 1]);
axis image
title('Coherent Change Detection (MLE)');
colormap gray
colorbar


%% Low-RCS/Shadow region identification
% Use entropy filtering
entropy_mask_im1 = entropyfilt(abs(im1),true(2*CCDwinSize+1));
entropy_mask_im2 = entropyfilt(abs(im2),true(2*CCDwinSize+1));
% the second argument is the neighborhood size to use

% Reconcile results from the two images
entropy_mask = min(entropy_mask_im1, entropy_mask_im2);

% Rescale to [0,1]
entropy_mask = rescale(entropy_mask);

% Plot the shadow regions
figure; 
imagesc(entropy_mask);
title('Entropy filtering (identify low-RCS shadow regions)')
axis image
colorbar


% Perform MLE change detection after masking shadow regions using this
% method
CCDimage = CCDbasic(im1,im2,CCDwinSize,'MLE',true);
% the last argument determines whether to implement low-RCS masking

% Display the CCD results
figure;
imagesc(CCDimage,[0 1]);
axis image
title('Coherent Change Detection (MLE; with entropy-based low-RCS masking)');
colormap gray
colorbar
