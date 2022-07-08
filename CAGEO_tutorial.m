%{
QDMlab: A MATLAB toolbox for analyzing quantum diamond microscope (QDM) 
magnetic field maps

Michael W.R. Volk (1), Roger R. Fu (1), Raisa Trubko (1,2), 
Pauli Kehayias (1,3), David R. Glenn (1,4) and Eduardo A. Lima (5)

(1) Department for Earth and Planetary Sciences, Harvard University, 20 Oxford Street, Cambridge, 02138, Massachusetts, USA
(2) Department of Physics, Worcester Polytechnic Institute, 100 Institute Road, Worcester, 01609, MA, USA
(3) Sandia National Laboratories, Albuquerque, 87123, NM, USA
(4) Department of Physics, Harvard University, 17 Oxford Street, Cambridge, 02138, Massachusetts, USA
(5) Department of Earth, Atmospheric, and Planetary Sciences, Massachusetts Institute of Technology, 77 Massachusetts Avenue - Bldg 54,Cambridge, 02139, Massachusetts, USA
%}

% A set of test data can be downloaded from 
% https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/MPMNZI 
% (Note: raw QDM data is rather big (several Gb). Downloads can take a long 
% time depending on the internet speed.)
close all; clc
%% Setting the data paths 
% Please set dataFolder to be the absolute path to the example data.
dataFolder = "C:\Users\micha\Dropbox\data\QDMlab_example_data";
nFolders = {fullfile(dataFolder, "data","complex_map_analysis", "MIL03346_AF100mT"),...
            fullfile(dataFolder, "data","complex_map_analysis", "MIL03346_AF28mT"),...
            fullfile(dataFolder, "data","complex_map_analysis","MIL03346_AF50mT"),...
            fullfile(dataFolder, "data","complex_map_analysis","MIL03346_NRM"),...
            fullfile(dataFolder, "data","dipole_inversion","Spherule_IRM1T"),...
            fullfile(dataFolder, "data","viscosity","ALH84001-FOV1_NOVISC_10-20-10"),...
            fullfile(dataFolder, "data","viscosity","ALH84001-FOV1_VISC_20-20")};

%% Measured images at different MW fields
% First lets have a look at what the images from the QDM camera look like
% at different microwave (MW) fields. This shows how the intensity of the
% fluorescence changes with aplied MW field.
% Lets load a file:
expData = load(fullfile(nFolders{4}, "run_00000.mat"));
% The for each of the 51 frequencies, we took a 1920x1200 image.
% Next the data has to be formatted from a 51x2304000 array to a 
% 1200x1920x51 array so we can plot it.

[binDataNorm, freq] = prepare_raw_data(expData, 1,1);
ax1 = subplot(1,2,1);
QDM_figure(binDataNorm(:,:,1), 'ax', ax1, 'led',1, 'title', sprintf('%.3f GHz',freq(1)))
ax2 = subplot(1,2,2);
QDM_figure(binDataNorm(:,:,25),'ax', ax2, 'led', 1, 'title', sprintf('%.3f GHz',freq(25)))

%% ODMR Raw Data Format
% Loading a run_*.mat file yields a structure with the following fields. 
% This measurement was done with 51 frequencies for the lower and upper 
% resonances each.
% NOTE: loading the data will take a little while

raw_data = load(fullfile(nFolders{1}, "run_00000.mat"))

%% Laser and reflected light measurements
% The QDM is capable of taking a reflected light image of the sample. This
% is extrem,ely useful for analizing the location of the magnetic sources.
% Also, a 'laser' image is usually saved. The 'laser' image only detects
% the fluorescence of the diamond, when there is no MW applied.
close all; clc;
fig = figure('Name', 'Laser and Refl. Light images', 'units', 'normalized', 'position', [0.2, 0.2, 0.7, 0.4]);
movegui(fig,'center')
data = load(fullfile(nFolders{1}, "4x4Binned", "final_fits_(4x4).mat"));
ax1 = subplot(1,2,1);
QDM_figure(data.ledImg, 'ax', ax1, 'led', true, 'title', 'reflected light');
ax2 = subplot(1,2,2);
QDM_figure(data.laserImg,'ax',ax2, 'title', 'laser', 'cbTitle','', 'unit', 'none');
linkaxes([ax1 ax2])

%% Mean Spectrum
% The mean measured (normalized) spectrum for this particular field 
% direction looks like:
close all; clc;
plot(raw_data.freqList, [raw_data.disp1, raw_data.disp2],'LineWidth',2)
title('Mean spectrum of all 2.3e6 pixel'); xlabel('f (Hz)'); ylabel('contrast')

%% Single Pixel Spectrum
%A single pixel on the other hand will be much more noisy and look 
% something like this:
close all; clc;

idx = randi(1000); %chose random pixel
plot(raw_data.freqList, [raw_data.imgStack1(:, idx); raw_data.imgStack2(:, idx)],...
    LineWidth=2)
title(['pixel #' num2str(idx)]); xlabel('f (Hz)'); ylabel('contrast')

% Each pixel needs to be fitted with either a triplet (N14) or doublet (N15) 
% on both sides of the spectrum. The fitting is done using ODMR_to_B111.

%% From ODMR spectra to B111 data
% Before you can calculate the B111 data from the ODMR spectra, you need to
% determine the global fluorescence fraction in the sample. 
globalFraction_estimator("binSize",1, "nRes", 1, "nRun", 1)

%% 
% This section shows you how to calculate the B_{111} data from the raw ODMR 
% spectra.
% NOTE: Running the following cell only works if you have a NVIDIA GPU 
%       installed in the computer. 
% NOTE: The execution of the following cell can take several minutes, to 
%       hours, depending on the CPU/GPU you are using. I suggest not 
%       executing it, since the data has already been calculated (i.e. 
%       see 1x1Binned/4x4Binned folders)

all_fits = ODMR_to_B111('nFolders', nFolders, ... % a cell of paths
     'binSizes', [1 4], ... % cell of integers. Here we chose 1 to show the unbinned data, and 4 (commonly used)
     'fieldPolarity', 0, ... % default: 0 = neg & pos; less common: 1 = neg, 2 = pos, (unusual 4 = nppn)
     'type', 2, ... % default: 2; fitting type, in most cases type = 2
     'globalFraction', 0.25, ... % default: 0.25; the ammount of global fluorescence as a ratio to the sample signal
     'checkPlot', false, ... % default: false; 
     'save', 1, ...
     'diamond', 'N14');

%% If no GPU is available
% We have pre-calculated the fits for all of the measurements. They are
% stored in the `1x1Binned` and `4x4Binned` folders within each measurement
% folder. The output of the line above will be similar to:

all_fits{1,1} = load(fullfile(dataFolder, "data","complex_map_analysis","MIL03346_NRM","1x1Binned", "final_fits_(1x1).mat"));
all_fits{1,2} = load(fullfile(dataFolder, "data","complex_map_analysis","MIL03346_NRM","4x4Binned", "final_fits_(4x4).mat"));
all_fits{2,1} = load(fullfile(dataFolder, "data","complex_map_analysis","MIL03346_AF28mT","1x1Binned", "final_fits_(1x1).mat"));
all_fits{2,2} = load(fullfile(dataFolder, "data","complex_map_analysis","MIL03346_AF28mT","4x4Binned", "final_fits_(4x4).mat"));
   
% If ODMR_to_B111 is called with n folders and m bins, ODMR_to_B111 will 
% return a (n x m) cell of the data that is ultimately stored in 
% final_fits_(mxm).mat.

%% Editing Maps
%% Subtract Blank
% Call subtract_blank without arguments. First, it will ask you to select 
% the measurement file (i.e. either B111dataToPlot.mat or 
% final_fits(nxn).mat). Second, you should select the blank file. 
% The blank is located in "dataFolder\data\blank_subtraction\Mar6_2020"

% subtract_blank();

% However, it is easier to call the function by passing nFiles and 
% blankFile keywords. nFiles i a cell of files, each of which will be 
% corrected for the diamond blank. The data is automatically saved 
% (`save`, true). A plot to check the diamond alignment (RECOMMENDED) 
% will be shown when checkPlot == true.

dataFile = fullfile(dataFolder, 'data', "viscosity","ALH84001-FOV1_NOVISC_10-20-10","4x4Binned","B111dataToPlot.mat");
blankFile = fullfile(dataFolder, 'data', "blank subtraction","Mar6_2020","4x4Binned","B111dataToPlot.mat");

subtract_blank( ...
    'nFiles', dataFile,...
    'blankFile', blankFile, ...
    'save', false, 'checkPlot', true);


%% Threshold Map
% To get rid of pixels that have not been calculated correctly 
% and show field values that are too high, threshold_map can be used. 
% Unfiltered, the map below shows pixels with field values > 0.1T, which 
% is outside of the mesurement region. Thus these pixels should be removed.

iFile = fullfile(dataFolder, "data","complex_map_analysis","MIL03346_NRM","4x4Binned","B111dataToPlot.mat");
dataToFilter = load(iFile);%%
QDM_figure(dataToFilter.B111ferro, 'preThreshold', 0, 'unit', 'microT');

%% Check Fits
% You can use `check_fits` to have a look at how well individual pixels
% were fitted.
% To do this, you need to load the final_fits_(nxn).mat and the 2 raw data
% files.

fits = load(fullfile(dataFolder, "data","complex_map_analysis","MIL03346_NRM","4x4Binned", "final_fits_(4x4).mat"));
raw1 = load(fullfile(dataFolder, "data","complex_map_analysis","MIL03346_NRM", "run_00000.mat"));
raw2 = load(fullfile(dataFolder, "data","complex_map_analysis","MIL03346_NRM", "run_00001.mat"));

% click on any point in the figure, to see the corresponding data/fits.
check_fits(fits, raw1, raw2)

%% threshold_QDM_map
% Using a threshold of 5G - meaning all pixels with B > 5G or B < -5G will 
% be replaced by the mean of the surrounding pixels - reveals the magnetic 
% signal of the rest of the sample. Thresholding large field values in 
% single pixels is especially useful, when the map needs to be converted 
% from B111 to Bz.

filteredData = threshold_QDM_map('nFiles', iFile, ...
    'threshold', 5);
QDM_figure(filteredData{1}.B111ferro, 'preThreshold', 0, 'unit', 'nT', ...
    'cbTitle','B_{111}');

%% Bz Conversion
% QDMdataprocessing automates the conversion process of B111 data into a 
% Bz map of remanent or induced signal (Fu et al., 2020).
QDMdataprocessing % unfiltered
%%
QDMdataprocessing % filtered
%%
data = load(fullfile(dataFolder, "data","complex_map_analysis", ...
                     "MIL03346_NRM","4x4Binned","Bz_uc0.mat"));
QDM_figure(data.Bz, 'unit','muT', 'scaleBar',500)

%% Analysis
%% Quantitative Analysis
% The QDM's ability to acquire reflected light images of the exact magnetic 
% field measurement region is  useful for locating ferromagnetic sources in 
% the sample. In the following example, a martian meteorite (MIL 03346, 
% Volk et al. 2020), the magnetization correcponds to the dark patches 
% (mesostasis) of the sample.

f = figure('Units','normalized', 'Position',[0,0,0.8,0.5])
movegui(f, 'center');

ax1 = subplot(1,2,1);
QDM_figure(filteredData{1}.B111ferro, 'preThreshold', 5, 'unit', 'T', ...
    'ax',ax1, 'title', 'Magnetic Data', 'xc', 1:1920, 'yc', 1:1200);
ax2 = subplot(1,2,2);
QDM_figure(filteredData{1}.ledImg, 'led', 1, 'unit', 'T', 'ax', ax2, ...
    'title', 'Reflected Light');
linkaxes([ax1 ax2])

%% Net Moment Analysis
%% Source Subtraction
% In some cases sources there are sources in our map that are not part of 
% the sample (dirt, dust). When we are trying to analyse individual sources 
% in the map these sources may have to be removed. As an example let us get 
% rid of the small source marked by the red circle.

fittingFile = fullfile(dataFolder, "data", "dipole_inversion", ...
    "Spherule_IRM1T","4x4Binned","Bz_uc0.mat");
fittingData = load(fittingFile);
QDM_figure(fittingData.Bz);
hold on
plot(368, 182, 'rO', 'MarkerSize',20)

%%
% subtract_source lets you select the top-left and bottom-right corners 
% of a rectangle that contains the source you want to remove.
sourceSubtracted = subtract_source('filePath', fittingFile, 'save', false);

QDM_figure(sourceSubtracted.Bz);
hold on
plot(266, 150, 'rO', 'MarkerSize',90)

%% Source Fitting
fit_sources_series(...
 fullfile(dataFolder, "data", "dipole_inversion", "Spherule_IRM1T","4x4Binned"), ...
 "upCont",{0,10,20});

%% complex_map_analysis
% Here we take a look at maps that do not show nice dipolar features like 
% the examples above. The signals are more complex (as in complicated not 
% the mathematical sense). The functions shown here allow us to extract 
% some information, even from these maps.

%% Demagnetization Behavior
% In paleomagnetism we often use stepwise thermal or alternating field
% demagnetization. We can estimate the 'thermal unblocking' or the
% coercivity of complex sources within the sample by analyzing changes in
% the magnetization patters. 
% With estimate_coercivity this can be done automatically.
% First create a cell with the path to the B111/Bz files.
nFolders = {fullfile(dataFolder, "data","complex_map_analysis", ...
                    "MIL03346_NRM", "4x4Binned");
            fullfile(dataFolder, "data","complex_map_analysis", ...
                    "MIL03346_AF28mT", "4x4Binned");
            fullfile(dataFolder, "data","complex_map_analysis", ...
                    "MIL03346_AF50mT", "4x4Binned");
            fullfile(dataFolder, "data","complex_map_analysis", ...
                    "MIL03346_AF100mT", "4x4Binned")};

% Then, you can execute the following line, and when prompted draw a 
% rectangle around the large complex source on the left. 
% Press esc to exit selection mode. 
% All images will get aligned with respect to the NRM so that small
% misorientations between measurements dont influence the analysis.
% A selection containing all pixels that exceed a field value of 0.1 *
% max(B)(NRM) and the number of pos. pixels within this selection is
% determined automatically.

demag_behavior(nFolders, 'selectionThreshold', 0.1, ...
    'fileName','B111dataToPlot.mat', 'checkPlot',1)

%% More than one source
% this can be done for multiple sources simultaneously by drawing
% rectangles around more than one source.

res = demag_behavior(nFolders, 'selectionThreshold', 0.1, ...
        'fileName','B111dataToPlot.mat', 'checkPlot',0);
demag_behavior_plot(res, "steps",[0, 28, 50, 100])

%% Viscous Magnetization
% Finally it is possible to visualize viscous parts of a sample. Here we
% have two measurements. One was measured with a viscous measurement
% protocol and a second one with a non-viscous (see full text for more info).
% The difference between the two gives the viscous magnetization.
% Similar to `estimate_coercivity`, the maps have to be aligned first. 
% The function `viscosity` does the alignement and then subtracts the two
% maps to show the viscous parts of the maps.

visc = viscosity("nonVisc", fullfile(dataFolder, "data", "viscosity", ...
                     "ALH84001-FOV1_NOVISC_10-20-10","4x4Binned", "Bz_uc0.mat"),...
                 "visc", fullfile(dataFolder, "data","viscosity", ...
                     "ALH84001-FOV1_VISC_20-20","4x4Binned", "Bz_uc0.mat"),...
                 "checkPlot", 1);
