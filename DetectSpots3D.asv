% [SpotList,Bg]=DetectSpots(InputImg,thresh,threshres,startsigma,method,ROISize,MaxSpotSize,MinInt)   : detects spots
% thresh : threshold in image to place new spots
% threshres : threshold in residuum to place new spots
% startsigma : starting sigma for fitting if one dimensional, than the
    % sigma defines the spots size. It can be either a number if the spot is
    % to be symmetrical in all dimensions or a vector of length same as the 
    % dimensionality of InputImg
% method : e.g 'idiv' method to use while fitting, see FitDataNDFastMask
% ROISize : X and Y radius of the ROI 
% MaxSpotSize : maximal size allowed for spots. It can be either a number
    % if the spot size limit is to be symmetrical in all dimensions or a 
    % vector of length same as the dimensionality of InputImg
% MinInt : minimu intensity allowed for spots


function [SpotList, Bg, fitted2, residual] = ...
    DetectSpots3D(InputImg, inputParam, Mask)

%% Find image dimensionality
numdim = ndims(InputImg);
numpar = 1 + numdim * 2;

%% Maximum number of iterrations allowed in each optimization round
if ~isfield(inputParam, 'maxiter')
    inputParam.maxiter = 500;
end

%% method : figure of merit to use for fitting 'mse', 'idiv' or 'fidiv'.
%     'mse' is least squares (assuming Gaussian noise).
%     'idiv' is Czesar's i-divergence (Sterling approximation).
%     'fidiv' is a fast version omitting data-dependent constants.
if ~isfield(inputParam, 'method')
    inputParam.method = 'idiv';
end

%% ROISize : X, Y(, Z) radius of the ROI 
if ~isfield(inputParam, 'ROISize')
    inputParam.ROISize = 18 * ones(1, numdim);
end

%% MaxSpotSize : maximal size allowed for spots. It can be either a number
%     if the spot size limit is to be symmetrical in all dimensions or a 
%     vector of length same as the dimensionality of InputImg
if ~isfield(inputParam, 'MaxSpotSize')
    inputParam.MaxSpotSize = 20;
end

%% MinInt : minimum intensity allowed for spots
if ~isfield(inputParam, 'MinInt')
    inputParam.MinInt = 10;
end

%% fixSigma: allows fixing spot size to increase the sensibility of the fit
if ~isfield(inputParam, 'fixSigma')
    inputParam.fixSigma = 0;
end

%% spotNrMax: maximum number of spots
if ~isfield(inputParam, 'spotNrMax')
    inputParam.spotNrMax = Inf;
end

%% maxVar: adaptive threshold coefficient selecting bright enough spots
if ~isfield(inputParam, 'maxVar')
    inputParam.maxVar = 50;
end

%% maxFits: Maximum number of DetectSpot cycles before it termintes to
%% prevent long time doing some nonsense
if ~isfield(inputParam, 'maxFits')
    inputParam.maxFits = 20;
end

%% Check if start sigma has the correct size
if length(inputParam.startsigma) ~= 1 && length(inputParam.startsigma) ~= numdim
    error('Parameter startsigma must be either a number or a vector of length equal to dimensionality of the input image');
end

if length(inputParam.startsigma) == 1
    inputParam.startsigma = repmat(inputParam.startsigma(1), 1, numdim);
end

%% Check if MaxSpotSize has the correct size
if length(inputParam.MaxSpotSize) ~= 1 && length(inputParam.MaxSpotSize) ~= numdim
    error('Parameter MaxSpotSize must be either a number or a vector of length equal to dimensionality of the input image');
end

if length(inputParam.MaxSpotSize) == 1
    inputParam.MaxSpotSize = repmat(inputParam.MaxSpotSize(1), 1, numdim);
end

%% Define gaussian filter sigma
switch numdim
    case 2
        filterS1 = max([6, 6; inputParam.startsigma]);
        filterS2 = max([4, 4; inputParam.startsigma]);
        filterS3 = max([2, 2; inputParam.startsigma]);
    case 3
        %filterS1 = max([6, 6, 3; inputParam.startsigma]);
        %filterS2 = max([4, 4, 2; inputParam.startsigma]);
        %filterS3 = max([2, 2, 1; inputParam.startsigma]);
        filterS1 = inputParam.startsigma;
        filterS2 = inputParam.startsigma / 3;
        filterS3 = inputParam.startsigma / 3;
end

%% Generate image mask
if nargin == 3
    useMask = gaussf(dip_image(Mask, 'double'), 3) * 2;
    useMask(useMask > 1) = 1;
end

%% Apply low frequency low pass filter
LLPfiltered = gaussf(InputImg, filterS1);
%filtered=gaussf(InputImg,max(6,startsigma));
%Bg = min(filtered);

%% Apply high frequency low pass filter
HLPfiltered = gaussf(InputImg, filterS2);

%% Apply bandpass filter
BPfiltered = HLPfiltered - LLPfiltered;

%filtered=gaussf(InputImg,max(4,startsigma));
%filteredorig = HLPfiltered;

%InputImg = BPfiltered;

Bg = mean(InputImg);

SpotList = [];
n = 1;
%BgFlag = 1;   % start with a fixed Background
BgFlag = 0;
oldSpotList = [];
activeImg = InputImg;
failedSpots = zeros(numdim, 0);

while 1
    %% Find local maxima with subpixel resolution
    % This function uses Gaussian interpolation to find the value and
    % position of the local maxima in the filtered image
    [Gravity, MaxVal] = dip_subpixelmaxima(BPfiltered * useMask, [], 'parabolic');
    
    %% Delete diverged maxima that are beyond the image boundary
    % The problem is that at the border the maxima tends to be fitted
    % as points somewhere outside of the boundix box ot the image these
    % need to be deleted
    in = find(sum(Gravity >= repmat(size(useMask), numel(MaxVal), 1) | Gravity <= 0, 2));
    Gravity(in, :) = [];
    MaxVal(in) = [];
    
    %% Delete spots which are outside the mask of the cell
    % The xyz coordinates are converted into a global index and this is
    % used to check if any of the maxima is where useMask is 0, those are
    % deleted
    % Do NOT use useMask here but the original mask. 
    in = find(logical(Mask(round(Gravity(:, 2)) + ...
        round(Gravity(:, 1)) * size(Mask, 2) + ...
        round(Gravity(:, 3)) * size(Mask, 1) * size(Mask, 2)) <= 0.5));
    Gravity(in, :) = [];
    MaxVal(in) = [];

    %% Sort the maxima
    % I will sort the maxima and 'plot' them and find those that are
    % outlying the straight line of random maxima around zero (those that
    % are substantially more bright than the others
    [spotInt, in] = sort(MaxVal'); %#ok<TRSRT>

    %% Fit a straight line through the intensities of 90% of the maxima
    %% leave out the weakest and strongest maxima
    x = floor(0.05 * numel(spotInt)) : ceil(0.95 * numel(spotInt));
    y = spotInt(x);
    [p, s, mu] = polyfit(x, y, 1);

    %% evaluate the linear prediction for all the points and the prediction
    %% of the standard deviation.
    x = 1 : numel(spotInt);
    [fy, d] = polyval(p, x, s, mu);
    %% Plot the linear model and the real intensity, the green line shows
    %% the 50% confidence interval
    if (0)
        semilogy(x,spotInt,'b.',x,fy,'r-',x,fy+numel(spotInt)/inputParam.maxVar*d,'g--', numel(spotInt)*0.05*[1,1], [1,1000],'k:', numel(spotInt)*0.95*[1,1], [1,1000],'k:',x,fy + numel(spotInt) / inputParam.maxVar * d,'m+')
        hl = legend('Local maxima', 'Linear regression', '50\times 50 % confidence interval', 'Central 90 % cutoff', 'selected cutoff', 'Location', 'SouthEast')
        xlabel('Local maximum index', 'FontSize', 20); ylabel('Local maximum value (log-scale)', 'FontSize', 20)
        % set(hl, 'FontSize', 18); set(gca, 'FontSize', 18); xlim([0,830])
        % print('/dev/shm/threshold.eps', '-depsc')
    end
    %% Select those points which lie well outside of the confidence
    %% interval, the 1/100 coefficient might need tweaking for different
    %% images
    in = fliplr(in(spotInt > fy + numel(spotInt) / inputParam.maxVar * d));
    Gravity = Gravity(in, :);
    MaxVal = MaxVal(in);
    
    %% Check if any of the spots are very similar to the failing spots
    % expand gravity matrix
    gm = repmat(Gravity', [1, 1, size(failedSpots, 2)]);
    % expand failed cell matrix
    fm = repmat(reshape(failedSpots, numdim, 1, size(failedSpots, 2)), [1, size(Gravity, 1), 1]);
    % compare them and find those which are less than 1 pixel from the
    % failed spot
    in = find(min(sqrt(sum((gm - fm) .^ 2, 1)), [], 3) < 1);
    % Remove spots which have failed fit last time
    Gravity(in, :) = [];
    MaxVal(in) = [];
    
    
    %% Scale the MaxVal to compensate for the intensity loss during
    %% filtering
    % The BPfiltered image has much lower intensity than the input image.
    % Therefore the spot intensity should be estimated from the real image.
    % It is done by a cubic interpolation of the input image at the center
    % of gravity acquired by the dip_subpixelmaxima earlier.

    if ~isempty(MaxVal)
        MaxVal = get_subpixel(activeImg, Gravity, 'cubic');
        % sometimes NaN value is given near the borded, hence these should
        % be discarded
        in = find(isnan(MaxVal));
        Gravity(in, :) = [];
        MaxVal(in) = [];
    end

    %% Terminate the loop if no more spots are found
    if isempty(MaxVal)
        fprintf('Found no more spots\n');
        
        %% If no spots were found at all, generate zero output variables
        if ~exist('fitted2', 'var')
            fitted2 = newim(size(InputImg));
        end
        if ~exist('residual', 'var')
            residual = InputImg;
        end
        break;
    end

    %% Check for the maximum number of spots allowed in the search
    % In the first search add only half of the maximum number of spots to
    % allow adding of more spots in the next rounds of fitting
    if n == 1
        % Generate an index which selects only the brightest spots until
        % the total number of spots is less than half the limit of number
        % of spots
        in = 1 : min(ceil(0.75 * inputParam.spotNrMax) - numel(SpotList) / numpar, numel(MaxVal));
    else
        % Generate an index which selects only the brightest spots until
        % the total number of spots is less than limit for number of spots
        in = 1 : min(inputParam.spotNrMax - numel(SpotList) / numpar, numel(MaxVal));
    end
    MaxVal = MaxVal(in)';
    Gravity = Gravity(in, :)';

    %% Generate a matrix with input parameters for all spots
    NewSpots = [MaxVal; Gravity; repmat(inputParam.startsigma', 1, numel(MaxVal))];
    
    %% variable for the number of fitted spots
    fprintf('Step %d, adding %d spots\n', n, size(NewSpots, 2));

    %% Form input vector with the start parameters of the fit
    initParms = [Bg, SpotList, NewSpots(:)'];
    
    %% Terminate the search if no spots are to be fitted
    if numel(initParms) < 1 + numpar
        fprintf('No spots found. Finishing search\n');
        % Generate output variables for no fit
        SpotList = [];
        Bg = 0;
        fitted2 = newim(size(InputImg, 1), size(InputImg, 2), size(InputImg, 3));
        residual = InputImg;
        break
    end

    %% Parameter mask
    % it is a vector of size 8:
    % 1:   background flag (0) variable background, (1) fixed background
    % 2:   spot maximum intensity (0) variable, (1) fixed, (2) global variable
    % 3-5: spot position (0) variable, (1) fixed, (2) global variable
    % 6-8: spot sigma (0) variable, (1) fixed, (2) global variable
    parmMask = [BgFlag, 0, zeros(1, numdim), inputParam.fixSigma * ones(1, numdim)];
    
    %% Perform the fit
    params = FitDataNDFastMask(initParms, parmMask, InputImg, numdim, ...
        inputParam.maxiter, inputParam.method, 0, inputParam.ROISize);
    
    %% Separate background offset and the spots parameters
    Bg = params(1);
    SpotList = params(2 : end);
    
    %% Check if any spot has deviated much from the starting point
    if any(initParms([3:7:end 4:7:end 5:7:end]) - params([3:7:end 4:7:end 5:7:end]) > 3)
        beep
    end

    %% reshape the spot parameters into a 7 x numSpots matrix
    SpotList = reshape(SpotList, [numpar, numel(SpotList) / numpar]);
    
    %% Fitted spot sigmas
	spotSigma = SpotList(2 + numdim : end, :);

    %% Find spots that are in any sigma exceeding the MaxSpotSize limit
    in1 = spotSigma >= repmat(size(useMask)', 1, size(SpotList, 2));

    %% Find spots that have any sigma zero or negative
    in2 = spotSigma <= 0;

    %% Fitted spot positions
    spotPos = SpotList(2 : 1 + numdim, :);

    %% Find spots that are outside of the field of view
    in3 = spotPos >= repmat(size(useMask)', 1, size(SpotList, 2)) | ...
        spotPos <= 0;
    
    %% Find spots that have maximum less than MinInt
    in4 = SpotList(1, :) < inputParam.MinInt;
    
    %% Combine all restrictions
    in = any(in1 | in2 | in3, 1) | in4;
    inI = find(~in);
    in = find(in);

    %% Find spots that are outside of the useMask and are not restricted
    %% otherwise
    if ~isempty(SpotList(:, inI))
        in1 = find(get_subpixel(useMask, SpotList(2 : numdim + 1, inI)', 'cubic') <= 0.5);
    else
        in1 = [];
    end
    
    %% Indices of unsatisfactory spots
    in = horzcat(in, inI(in1)); %#ok<AGROW>
    if any(in)
        beep
    end

    %% Remember unsatisfactory spots to avoid them next time
    failedSpots = horzcat(failedSpots, ...
        NewSpots(2 : 1 + numdim, ...
        in(in > size(SpotList, 2) - size(NewSpots, 2)) - ...
        (size(SpotList, 2) - size(NewSpots, 2)))); %#ok<AGROW>
    
    %% Remove all unsatisfactory spots
    SpotList(:, in) = [];
    fprintf('Deleted %d spots\n', numel(in) + numel(in1))
    
    %% Reshape the SpotList into a vector of input parameters for fit
    SpotList = SpotList(:)';
    initParms = [Bg, SpotList];

    %% Terminate the search if no spots are to be fitted
    if numel(initParms) < 1 + numpar
        fprintf('No spots found. Finishing search\n');
        SpotList = [];
        Bg = 0;
        fitted2 = newim(size(InputImg, 1), size(InputImg, 2), size(InputImg, 3));
        residual = InputImg;
        break
    end

    %% From the second run on, the background is fitted as well
    % althought I fit the background from the beginning anyway
    BgFlag = 0;

    %% Peform the fit
    [params, junk, fitted2, residual] = ...
        FitDataNDFastMask(initParms, parmMask, InputImg, numdim, ...
        inputParam.maxiter, inputParam.method, 0, inputParam.ROISize);

    %% Separate the bacground value from the spot parameters
    Bg = params(1);
    SpotList = params(2 : end);

    %% The program can lock itself in infinite loop when a spot is deleted
    % for not fulfilling a certain condition and then this spot is added
    % again in the next spot search. Therefore we need to compare
    % successive search results and terminate if they don't differ from
    % each other    
    if numel(oldSpotList) == numel(SpotList)
        if sum(abs(oldSpotList - SpotList)) < 1
            fprintf('Found no more spots, removed unfulfilling spots\n');
            break;
        end
    end
    oldSpotList = SpotList;
    residual = InputImg - fitted2; % To make it work for ROIs as well

    %% Terminate run if the maximum number of spots has been reached
    if numel(SpotList) / numpar == inputParam.spotNrMax
        fprintf('Reached the maximum allowed number of spots\n')
        break
    end

    %% Terminate if the number of iterations is exceeded
    if n >= inputParam.maxFits
        fprintf('Reached the maximum number of attempted fits\n')
        break
    end

    %% Increment the number of fits counter
    n = n + 1;

    %% Fit the residual in the next round
    activeImg = residual;
    
    %% Filter with low frequency low pass filter
    LLPfiltered = gaussf(residual, filterS1);
    %% Filter with high frequency low pass filter
    HLPfiltered = gaussf(residual, filterS3);

    %% Apply bandpass filter
    BPfiltered = HLPfiltered - LLPfiltered;
end
SpotList = reshape(SpotList, numpar, length(SpotList) / numpar);
SpotList = SpotList';
[junk, i] = sort(SpotList(:, 1), 'descend');
SpotList = SpotList(i, :);

fprintf('Number of spots found: %d\n', size(SpotList, 1))

% tmp=SpotList(2,:);   % Swap X and Y coordinates of position
% SpotList(2,:)=SpotList(3,:);
% SpotList(3,:)=tmp;
% 
% tmp=SpotList(4,:);   % Swap X and Y coordinates of sizes
% SpotList(4,:)=SpotList(5,:);
% SpotList(5,:)=tmp;
