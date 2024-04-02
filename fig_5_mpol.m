close all; clear; clc

% Parts of this code are based on the Matlab script files included as part 
% of the GOTCHA CCD challenge problem dataset
%
% The original version of the scripts were written by 
% Steven Scarborough and LeRoy Gorham (AFRL/RYAP)
% Email:  steven.scarborough@wpafb.af.mil / leroy.gorham@wpafb.af.mil

% for pretty pictures
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex');


% NOTE: This script has a slightly longer runtime; expect 1-2 minutes 
% before seeing results 


%% Problem Parameters

% Define the path to the base directory of the dataset
datadir = '~/Downloads/ChangeDetectionDataset_GOTCHA/CCD-CP-XPol';

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


% Apply normalization?
applyNormalization = true;


%% Read in data 
% Determine the file names of the input files
% ... for the HH polarization
im1filename = sprintf('%s/HH/%s/c00007a283p50.mat',datadir,pass1);
im2filename = sprintf('%s/HH/%s/c00007a283p50.mat',datadir,pass2);

% Load in the reference and missions images
data1 = load(im1filename);
data2 = load(im2filename);

% Crop the image to requested size
im1_HH = data1.SARdataOut(minRow:maxRow,minCol:maxCol);
im2_HH = data2.SARdataOut(minRow:maxRow,minCol:maxCol);


% ... now repeat for the VV polarization
im1filename = sprintf('%s/VV/%s/c00007a283p50.mat',datadir,pass1);
im2filename = sprintf('%s/VV/%s/c00007a283p50.mat',datadir,pass2);

% Load in the reference and missions images
data1 = load(im1filename);
data2 = load(im2filename);

% Crop the image to requested size
im1_VV = data1.SARdataOut(minRow:maxRow,minCol:maxCol);
im2_VV = data2.SARdataOut(minRow:maxRow,minCol:maxCol);


% ... and for the HV polarization
im1filename = sprintf('%s/HV/%s/c00007a283p50.mat',datadir,pass1);
im2filename = sprintf('%s/HV/%s/c00007a283p50.mat',datadir,pass2);

% Load in the reference and missions images
data1 = load(im1filename);
data2 = load(im2filename);

% Crop the image to requested size
im1_HV = data1.SARdataOut(minRow:maxRow,minCol:maxCol);
im2_HV = data2.SARdataOut(minRow:maxRow,minCol:maxCol);


%% Normalization of the different polarization images

% % To test effectiveness of the normalization procedure
% im1_HH = 10*im1_HH;
% im1_VV = 10*im1_VV;
% im1_HV = 10*im1_HV;


if(applyNormalization)    
    % Now, we normalize across images
    % Use the Frobenius norm
    norm_im1_HH = norm(im1_HH, 'fro');
    norm_im1_VV = norm(im1_VV, 'fro');
    norm_im1_HV = norm(im1_HV, 'fro');

    norm_im2_HH = norm(im2_HH, 'fro');
    norm_im2_VV = norm(im2_VV, 'fro');
    norm_im2_HV = norm(im2_HV, 'fro');
    
    % Scale all images to have the same norm as the time-1, HH image
    im1_VV = (norm_im1_HH/norm_im1_VV)*im1_VV;
    im1_HV = (norm_im1_HH/norm_im1_HV)*im1_HV;

    im2_HH = (norm_im1_HH/norm_im2_HH)*im2_HH;
    im2_VV = (norm_im1_HH/norm_im2_VV)*im2_VV;
    im2_HV = (norm_im1_HH/norm_im2_HV)*im2_HV;
end


%% Mean subtraction
nbr_sz = (2*CCDwinSize+1)^2;        % Neighborhood size
im1_HH = im1_HH - conv2(im1_HH, ones(2*CCDwinSize+1), 'same')/nbr_sz;
im1_VV = im1_VV - conv2(im1_VV, ones(2*CCDwinSize+1), 'same')/nbr_sz;
im1_HV = im1_HV - conv2(im1_HV, ones(2*CCDwinSize+1), 'same')/nbr_sz;

im2_HH = im2_HH - conv2(im2_HH, ones(2*CCDwinSize+1), 'same')/nbr_sz;
im2_VV = im2_VV - conv2(im2_VV, ones(2*CCDwinSize+1), 'same')/nbr_sz;
im2_HV = im2_HV - conv2(im2_HV, ones(2*CCDwinSize+1), 'same')/nbr_sz;


%% Compute change map using the multi-polarization test statistic 
% (without low-RCS masking)
CCDimage = CCD_mpol( ...
    cat(3,im1_HH,im1_VV,im1_HV), ...    % all polarizations for time-1 signal 
    cat(3,im2_HH,im2_VV,im2_HV), ...    % all polarizations for time-1 signal 
    CCDwinSize ...                      % neighborhood size
                );


% Display the CCD results
figure
imagesc(CCDimage,[0 1]);
axis image
title('Coherent Change Detection (Multipolar; HH/VV/HV polarizations)');
colormap gray
colorbar