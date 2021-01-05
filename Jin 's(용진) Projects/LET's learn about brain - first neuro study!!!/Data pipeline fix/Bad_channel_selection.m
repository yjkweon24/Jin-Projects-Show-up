%{
Qs

1. 

2. 

%}


%% Happe

%{
Need to make specific inputs for channels as well:
    https://github.com/lcnhappe/happe/blob/master/HAPPE_v1.0_README.pdf
 
 Script:
    https://github.com/lcnhappe/happe/blob/master/HAPPE_pipeline_v1_0.m
%}


% 3. list channels of interest, including the 10-20 channels. User defined channels occur at the end of the sequence e.g. 'E39' 'E40'
%the 18 "10-20" channels that NEED to be in the chan_IDs: 'FP1' 'FP2' 'F3'
% 'F4' 'F7' 'F8' 'C3' 'C4' 'T3' 'T4' 'PZ' 'O1' 'O2' 'T5' 'T6' 'P3' 'P4' 'Fz'
chan_IDs={'FP1' 'FP2' 'F3' 'F4' 'F7' 'F8' 'C3' 'C4' 'T3' 'T4' 'PZ' 'O1' 'O2' ...
    'T5' 'T6' 'P3' 'P4' 'Fz' 'E27' 'E23' 'E19' 'E20' 'E28' 'E13' 'E41' 'E40' 'E46'...
    'E47' 'E75' 'E3' 'E4' 'E123' 'E118' 'E112' 'E117' 'E109' 'E102' 'E98' 'E103'};

chan_index=[1:length(chan_IDs)];



% do you want to do segment rejection using all user-specified channels above ( = 0) or a subset of channels in an ROI ( = 1)?
ROI_channels_only = 0;


% if you want to do ROI segment rejection, which channels should be used in the ROI?
% ex ROI_channels = {'E27','E20'};
ROI_channels = {'E27','E20'};
% if you are referencing to another channel/subset of channels, what are they?
% make sure to use the channel name given above
% ex ROI_channels = {'E57','E100'};



% Normed joint probability of the average log power from 1 to 125 Hz across channels. Z > 3 are removed.  Performed twice
%pop_rejchan & eeg_checkset are functions from EEGLAB.

%% crude bad channel detection using spectrum criteria and 3SDeviations as channel outlier threshold, done twice
% all the functions:
% https://github.com/lcnhappe/happe/tree/master/Packages/eeglab14_0_0b/functions/popfunc 
% https://github.com/lcnhappe/happe/tree/master/Packages/eeglab14_0_0b/functions/adminfunc

    EEG = pop_rejchan(EEG, 'elec',chan_index,'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 125]);
    EEG.setname='rawEEG_f_cs_ln_badc';
    EEG = eeg_checkset( EEG );
    
    EEG = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 125]);
    EEG.setname='rawEEG_f_cs_ln_badc2x';
    EEG = eeg_checkset( EEG );
    selected_channel_locations=EEG.chanlocs;
    
    %save the names of the rejected channels for output table after the pipeline finishes
    selected_channel_labels={selected_channel_locations.labels};
    bad_channels_removed= setdiff(chan_IDs, selected_channel_labels);
    [~,ROI_indices_in_selected_chanlocs] = intersect(selected_channel_labels,ROI_channels);
    
%% interpolate the channels that were flagged as bad earlier:
    EEG = pop_interp(EEG, full_selected_channels, 'spherical');
    EEG.setname='wavcleanedEEG_ICA_MARA_rej_chan_int';
    EEG = eeg_checkset(EEG );     
     
     
     
     
     
    
    
    
%% PREP 1, 2, 5 

%{
https://github.com/VisLab/EEG-Clean-Tools/blob/master/PrepPipeline/utilities/findNoisyChannels.m
    
https://github.com/VisLab/EEG-Clean-Tools/blob/master/PrepPipeline/utilities/findUnusableChannels.m

https://github.com/VisLab/EEG-Clean-Tools/blob/master/PrepPipeline/prepPostProcess.m

https://github.com/VisLab/EEG-Clean-Tools/blob/master/PrepPipeline/utilities/interpolateChannels.m
    
https://github.com/VisLab/EEG-Clean-Tools/blob/master/PrepPipeline/utilities/private/spherical_interpolate.m
%}



%% Method 1: Unusually high or low amplitude (using robust std) - DEviation criteria 
channelDeviation = 0.7413 *iqr(data); % Robust estimate of SD
channelDeviationSD =  0.7413 * iqr(channelDeviation);
channelDeviationMedian = nanmedian(channelDeviation);
noisyOut.robustChannelDeviation(evaluationChannels) = ...
    (channelDeviation - channelDeviationMedian) / channelDeviationSD;

% Find channels with unusually high deviation 
badChannelsFromDeviation = ...
    abs(noisyOut.robustChannelDeviation) > ...
             noisyOut.robustDeviationThreshold | ...
             isnan(noisyOut.robustChannelDeviation);
badChannelsFromDeviation = originalChannels(badChannelsFromDeviation);
noisyOut.noisyChannels.badChannelsFromDeviation = badChannelsFromDeviation(:)';
noisyOut.channelDeviationMedian = channelDeviationMedian;
noisyOut.channelDeviationSD = channelDeviationSD;





%% Method 2: Global correlation criteria (from Nima Bigdely-Shamlo)
channelCorrelations = ones(WCorrelation, numberChannels);
noiseLevels = zeros(WCorrelation, numberChannels);
channelDeviations = zeros(WCorrelation, numberChannels);
n = length(correlationWindow);
xWin = reshape(X(1:n*WCorrelation, :)', numberChannels, n, WCorrelation);
dataWin = reshape(data(1:n*WCorrelation, :)', numberChannels, n, WCorrelation);
parfor k = 1:WCorrelation 
    eegPortion = squeeze(xWin(:, :, k))';
    dataPortion = squeeze(dataWin(:, :, k))';
    windowCorrelation = corrcoef(eegPortion);
    abs_corr = abs(windowCorrelation - diag(diag(windowCorrelation)));
    channelCorrelations(k, :)  = quantile(abs_corr, 0.98);
    noiseLevels(k, :) = mad(dataPortion - eegPortion, 1)./mad(eegPortion, 1);
    channelDeviations(k, :) =  0.7413 *iqr(dataPortion);
end;
dropOuts = isnan(channelCorrelations) | isnan(noiseLevels);
channelCorrelations(dropOuts) = 0.0;
noiseLevels(dropOuts) = 0.0;
clear xWin;
clear dataWin;
noisyOut.maximumCorrelations(evaluationChannels, :) = channelCorrelations';
noisyOut.noiseLevels(evaluationChannels, :) = noiseLevels';
noisyOut.channelDeviations(evaluationChannels, :) = channelDeviations';
noisyOut.dropOuts(evaluationChannels, :) = dropOuts';
thresholdedCorrelations = ...
    noisyOut.maximumCorrelations < noisyOut.correlationThreshold;
fractionBadCorrelationWindows = mean(thresholdedCorrelations, 2);
fractionBadDropOutWindows = mean(noisyOut.dropOuts, 2);

% Remap channels to their original numbers
badChannelsFromCorrelation = find(fractionBadCorrelationWindows > noisyOut.badTimeThreshold);
noisyOut.noisyChannels.badChannelsFromCorrelation = badChannelsFromCorrelation(:)';
badChannelsFromDropOuts = find(fractionBadDropOutWindows > noisyOut.badTimeThreshold);
noisyOut.noisyChannels.badChannelsFromDropOuts = badChannelsFromDropOuts(:)';
noisyOut.medianMaxCorrelation =  median(noisyOut.maximumCorrelations, 2);

% Bad so far by amplitude and correlation (take these out before doing ransac)
noisyChannels = union(noisyOut.noisyChannels.badChannelsFromDeviation, ...
    union(noisyOut.noisyChannels.badChannelsFromCorrelation, ...
          noisyOut.noisyChannels.badChannelsFromDropOuts));





%% Detect constant or NaN channels and remove from consideration
nanChannelMask = sum(isnan(data), 1) > 0;
noSignalChannelMask = mad(data, 1, 1) < 10e-10 | std(data, 1, 1) < 10e-10;
noisyOut.noisyChannels.badChannelsFromNaNs = evaluationChannels(nanChannelMask);
noisyOut.noisyChannels.badChannelsFromNoData = evaluationChannels(noSignalChannelMask);
evaluationChannels = setdiff(evaluationChannels, ...
    union(noisyOut.noisyChannels.badChannelsFromNaNs, ...
          noisyOut.noisyChannels.badChannelsFromNoData));
data = signal.data;
data = double(data(evaluationChannels, :))';  
[signalSize, numberChannels] = size(data);



%% Detect channels that are constant for large periods or have NaNs
function [nanChannels, noSignalChannels] = ...
                        findUnusableChannels(signal, evaluationChannels) 

nanChannels = sum(isnan(signal.data(evaluationChannels, :)), 2) > 0;
noSignalChannels = mad(signal.data(evaluationChannels, :), 1, 2) < 10e-10 | ...
    std(signal.data(evaluationChannels, :), 1, 2) < 10e-10;

% Now convert to channel numbers
nanChannels = evaluationChannels(nanChannels);
noSignalChannels = evaluationChannels(noSignalChannels);



%% Combine bad channels detected from all methods
noisy = noisyOut.noisyChannels;
noisyOut.noisyChannels.badChannelsFromLowSNR = ...
    intersect(noisy.badChannelsFromHFNoise, noisy.badChannelsFromCorrelation);
noisyChannels = union(noisyChannels, ...
    union(union(noisy.badChannelsFromRansac, ...
          noisy.badChannelsFromHFNoise), ...
    union(noisy.badChannelsFromNaNs, ...
          noisy.badChannelsFromNoData)));
noisyOut.noisyChannels.all = noisyChannels(:)';
noisyOut.medianMaxCorrelation =  median(noisyOut.maximumCorrelations, 2);
 






%% Interpolation

function signal = interpolateChannels(signal, targetChannels, sourceChannels)
% Interpolate the targetChannels rows of signal.data using spherical splines
%
% signal = interpolateChannels(signal, targetChannels)
% [signal, W] = interpolateChannels(signal, targetChannels, sourceChannels)
%
% Input:
%   signal      Structure with data and chanlocs fields compatible with
%               an EEGLAB EEG structure (requires .data and .chanlocs fields)
%   targetChannels  Vector of channel numbers for channels to interpolate (required)
%   sourceChannels  Vector of channel numbers to use in the interpolation
%               (default all channels except destination channels)
%
% Output:
%   signal      Input signal structure with the dest_chans rows of 
%               signal.data replaced with their interpolated values
%
%
% Written by Kay Robbins 8/24/2014 UTSA 
%

% Check the parameters
if nargin < 1 || ~isstruct(signal) || ~isfield(signal, 'data') || ...
   ~isfield(signal, 'chanlocs')      
    error('interpolateChannels:NotEnoughArguments', ...
          'first argument must be a structure with data and chanlocs fields');
elseif ~exist('targetChannels', 'var') || isempty(targetChannels) || ...
    min(targetChannels) < 0 || max(targetChannels) > size(signal.data, 1)  
    error('interpolateChannels:InterpolatedChannelsOutOfRange', ...
          'interpolated (dest) channels must be in dim 1 of signal.data');
elseif ~exist('sourceChannels', 'var') || isempty(sourceChannels)
    sourceChannels = 1:size(signal.data, 1); 
end
sourceChannels = setdiff(sourceChannels, targetChannels);  % Src and dest not same
if isempty(sourceChannels) || min(sourceChannels) < 0 || max(sourceChannels) > size(signal.data, 1)
    error('interpolateChannels:SourceChannelsOutOfRange', ...
          'source channels must be in dim 1 of signal.data');
end

% Perform the interpolation -- currently using EEGLAB eeg_interp
sourceChannels = sort(sourceChannels);
targetChannels = sort(targetChannels);
channelsConsidered = union(sourceChannels, targetChannels);
[~, reindexedTargetChannels] = intersect(channelsConsidered, targetChannels);
signaltemp = signal;
signaltemp.data = signal.data(channelsConsidered, :);
signaltemp.chanlocs = signal.chanlocs(channelsConsidered);
signaltemp.nbchan = length(signaltemp.chanlocs);
signaltemp = eeg_interp(signaltemp, reindexedTargetChannels, 'spherical');
reindex = false(1, length(signal.chanlocs));
reindex(targetChannels) = true;
reindex = reindex(channelsConsidered);
signal.data(targetChannels, :) = signaltemp.data(reindex, :);




%% Spherical interpolation

function [W, Gss, Gds, Hds] = ...
     spherical_interpolate(src, dest, lambda, order, type, tol)
% Caclulate an interpolation matrix for spherical interpolation
%
% W = sphericalSplineInterpolate(src, dest, lambda, order, type, tol)
%
% Inputs:
%  src    - [3 x N] old electrode positions
%  dest   - [3 x M] new electrode positions
%  lambda - [float] regularisation parameter for smoothing the estimates (1e-5)
%  order  - [float] order of the polynomial interpolation to use (4)
%  type - [str] one of; ('spline')
%             'spline' - spherical Spline 
%             'slap'   - surface Laplacian (aka. CSD)
%  tol    - [float] tolerance for the legendre poly approx (1e-7)
%
% Outputs:
%  W      - [M x N] linear mapping matrix between old and new co-ords
%
% Based upon the paper: Perrin89

% Copyright 2009-     by Jason D.R. Farquhar (jdrf@zepler.org)
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents, provided this copyright
% notice is retained, and note is made of any changes that have been
% made. This software and documents are distributed without any
% warranty, express or implied. 
%
% Modified by Kay Robbins 8/24/2014: Minor cleanup and simplification
% Warning --- still in progress

if ( nargin < 3 || isempty(lambda) ) lambda = 1e-5; end; %#ok<SEPEX>
if ( nargin < 4 || isempty(order) ) order = 4; end; %#ok<SEPEX>
if ( nargin < 5 || isempty(type)) type = 'spline'; end; %#ok<SEPEX>
if ( nargin < 6 || isempty(tol) ) tol = eps; end; %#ok<SEPEX>

% Map the positions onto the sphere (not using repop, by JMH)
src  = src./repmat(sqrt(sum(src.^2)), size(src, 1), 1);  
dest = dest./repmat(sqrt(sum(dest.^2)), size(dest, 1), 1); 

% Calculate the cosine of the angle between the new and old electrodes. If
% the vectors are on top of each other, the result is 1, if they are
% pointing the other way, the result is -1
cosSS = src'*src;  % angles between source positions
cosDS = dest'*src; % angles between destination positions

% Compute the interpolation matrix to tolerance tol
[Gss] = interpMx(cosSS, order, tol);       % [nSrc x nSrc]
[Gds, Hds] = interpMx(cosDS, order, tol);  % [nDest x nSrc]

% Include the regularisation
if lambda > 0 
    Gss = Gss + lambda*eye(size(Gss));
end

% Compute the mapping to the polynomial coefficients space % [nSrc+1 x nSrc+1]
% N.B. this can be numerically unstable so use the PINV to solve..
muGss = 1;  % Used to improve condition number when inverting. Probably unnecessary
C = [Gss  muGss*ones(size(Gss, 1),1); muGss*ones(1, size(Gss,2)) 0];
iC = pinv(C);

% Compute the mapping from source measurements and positions to destination positions
if ( strcmpi(type, 'spline') )
  W = [Gds ones(size(Gds, 1), 1).*muGss]*iC(:, 1:end-1); % [nDest x nSrc]
elseif (strcmpi(type, 'slap'))
  W = Hds*iC(1:end-1, 1:end-1); % [nDest x nSrc]
end
return;
%--------------------------------------------------------------------------
function [G, H]=interpMx(cosEE, order, tol)
% compute the interpolation matrix for this set of point pairs
if ( nargin < 3 || isempty(tol) ) 
    tol = 1e-10; 
end
G = zeros(size(cosEE)); H = zeros(size(cosEE));
for i = 1:numel(cosEE);
   x = cosEE(i);
   n = 1; Pns1 = 1; Pn = x;       % seeds for the legendre ploy recurrence
   tmp  = ( (2*n + 1) * Pn ) / ((n*n + n).^order);
   G(i) = tmp ;         % 1st element in the sum
   H(i) = (n*n + n)*tmp;  % 1st element in the sum
   dG = abs(G(i)); dH = abs(H(i));
   for n = 2:500; % do the sum
      Pns2 = Pns1; Pns1 = Pn; 
      Pn=((2*n - 1)*x*Pns1 - (n - 1)*Pns2)./n; % legendre poly recurrence
      oGi = G(i);  oHi = H(i);
      tmp = ((2*n+1) * Pn) / ((n*n+n).^order);
      G(i) = G(i) + tmp;             % update function estimate, spline interp     
      H(i) = H(i) + (n*n + n)*tmp;   % update function estimate, SLAP
      dG = (abs(oGi - G(i)) + dG)/2; 
      dH = (abs(oHi - H(i)) + dH)/2; % moving ave gradient est for convergence
      if (dG < tol && dH < tol) 
          break; 
      end;           % stop when tol reached
   end
end
G = G./(4*pi);
H = H./(4*pi);
return;



%% Remove interpolated channels if requested
 try
    interpolatedChannels = EEG.etc.noiseDetection.reference.interpolatedChannels.all;
    if postOut.removeInterpolatedChannels && ~isempty(interpolatedChannels)
      EEG.etc.noiseDetection.postProcess.removedChannelNumbers = ...
           interpolatedChannels; 
      EEG.etc.noiseDetection.postProcess.removedChannelData = ...
           EEG.data(interpolatedChannels, :);
      EEG.etc.noiseDetection.postProcess.removedChanlocs = ...
           EEG.chanlocs(interpolatedChannels);
      EEG.chanlocs(interpolatedChannels) = [];
      EEG.data(interpolatedChannels, :) = [];
      EEG.nbchan = length(EEG.chanlocs);                              
    else
       EEG.etc.noiseDetection.postProcess.removedChannelNumbers = [];
       EEG.etc.noiseDetection.postProcess.removedChannelLocs = [];
    end   
catch mex
    errorMessages.removeInterpolated = ...
        ['postProcessing failed to remove interpolated channels: ' ...
          getReport(mex, 'basic', 'hyperlinks', 'off')];
    EEG.etc.noiseDetection.errors.status = 'post-processing error';
    EEG.etc.noiseDetection.errors.postProcess = errorMessages;
    if strcmpi(postOut.errorMsgs, 'verbose')
        warning('[%s]\n%s', mex.identifier, ...
            getReport(mex, 'extended', 'hyperlinks', 'on'));
    end
    return;
 end
 postOut = EEG.etc.noiseDetection.postProcess;
