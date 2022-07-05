%% SIMPLE STEP-BY-STEP INSTRUCTIONS FOR BASIC FITTING
% This is a short introduction for fitting and calculating a Bz map fdrom
% raw ODMR data. For more information check the function help and the paper
%
% 
%% clear all Matlab variables and plots beforehand
% This is only for convenience, not needed.
clear all; close all; clc

%% Data location
% Set the folder location for the tutorial data. As an excample: my data 
% is stored at:
% "C:\Users\micha\OneDrive\Desktop\CAGEO_data"
% You need to change this to the location you downloaded the data to.

dataFolder = "C:\Users\micha\OneDrive\Desktop\CAGEO_data";

%% Check Data
% Here we have a quick look at what the data looks like. 
% The way QDMlab specifies the ODMR data is to separate the full spectrum 
% into a 'high' and 'low' frequency part. 
% check_ODMR allows you to check how the spectra look like at certain
% positions of the image. Just click on either of the top two (reflected
% light or laser) images and the spectra of the closest pixel will be
% shown.

check_ODMR(dataFolder)

%% determine the globalFraction
% First we need to load one of the two run_000*.mat files. They contain the
% data for the two field polarities. Either one is fine in this case.

expData = load(fullfile(dataFolder, 'run_00000.mat'));

% Using this data we can have a look at how much the global Illumination
% plays a role in this dataset. The example data provided in the supplement
% is already binned by 4. This means the 'hot pixels' have mostly been
% averaged out, already and it is more difficult to see how much of the 
% global signal needs to be corrected.
% Try using the bottom slider, from experience, a value around 0.35 seems 
% to be a good value.

globalFraction_estimator(expData, "binSize", 1, 'nRes', 1, 'nRun', 1);
globalFraction = 0.35;

%% Calculating the B111 fields from the raw spectra 
% The high and low frequency raw data for each pixel and field direction 
% needs to be fitted. This is done with ODMR_to_B111. 
% The sample is very magnetic, so that a number of pixels will not be
% fitted correctly. The ODMR peaks in some cases lie outside of the
% frequency range of the measurement. 
% If 'checkPlot' is set to 1, you can click on the high X2 pixels to
% check the fits for the pixels that are not fitted well.
% NOTE: Running this will only work with a working NVIDIA GPU.

fits = ODMR_to_B111("nFolders", dataFolder, 'globalFraction', globalFraction, ...
    "binSizes", 1, "checkPlot", 1, "diamond", 'N14');

%% Alternatively, if you can not run this yourself, the
% 'final_fits_(1x1).mat' is also in the '1x1Binned' folder, where you
% stored the data. Execute the following line to check its contents.
fits = {load(fullfile(dataFolder, '1x1Binned', 'final_fits_(1x1).mat'))};

%% Fits structure 
% Finally the data is aved to the disk, but can also be returned. Here we
% use the variable 'fits'. This contains the same informationa s the
% final_fits_(nxn).mat files. You can pass a more than one folder into
% ODMR_to_B111, therefore the variable 'fits' is a cell containing all (if
% more than one) fits. Here we only have one, so the first element of 'fits' 
% is the correct one.
% It contains all the information needed for further analysis (e.g. LED 
% and laser images, pixelAlerts, B111 field data).

fits{1}

%% Removing hot pixels for Bz conversion
% As seen earlier, the data here contains some pixels, that have not been
% fitted properly, because the peaks lie outside of the frequency range of 
% the measurement. Therefore, we need to remove these hot pixels before
% transforming the B111 data to the more convenient Bz.
iFile = fullfile(dataFolder, '1x1Binned', 'final_fits_(1x1).mat');
% We are replacing all pixels with field values > 3G (300 µT). By changing the
% 'treshold' value you can remove more or less hot pixels.
threshold_QDM_map('nFiles', iFile, 'threshold', 3, 'save', true, 'checkPlot', 1);

%% converting to Bz
% To convert the B111 data into Bz, 'QDMdataprocessing' is used. Select the
% '_thresh' map we saved in the previous step. 
% This function will ask several things:
%   - Do you want to calculate the para- (induced) or 
%        ferro-magnetic (remanent) component. 
%   - Pixel size:
%       which has to be determined from you mircoscopes optical 
%       path. The test data provided has a pixel size of 4.7 µm.
%   - NV-sample distance: Estimate of the sample to sensor distance
% Then it opens a plot where you can crop the data to your liking. Press
% 'S' to skip cropping.
% After cropping, the Bz plot opens and lets you do some basic background
% subtraction. By clicking on a 'background' pixel, the field value of that
% pixel will be subtracted from the whole image.
% Finally there is the option for a more sophisticated background
% subtraction, interpolation, and upward continuation.
%
% if you are interested you can perform the sdame conversion on the
% 'unfiltered' data and have a look at the difference.

QDMdataprocessing
