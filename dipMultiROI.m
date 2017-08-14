% function [ROIs,points,n]=dipMultiROI(img, numClasses) : Allows the user to select multiple ROIs
% img : Input image for ROI selection
% numClasses : if not given or empty all ROIs will automatically correspond to new labels. If numClasses is given, all consecutive ROIs will belong to the same label until empty ROI is detected when the next label starts.
%
% ROIs : label image with ROIs numbered
% points : cell array of clicked locations
% n : number of ROIs selected by user
%
% Example:
% img=readim;
% [ROIs,points,n]=dipMultiROI(img);
% m=measure(ROIs,img,'Sum');
% m

function [ROIs,points,n]=dipMultiROI(img, numClasses)
if nargin < 2
    numClasses=[];
end
if ndims(img) > 2
    error('dipMultiROI only works for 2D images. Try "squeeze" if your image is already 2D');
end
ROIs=newim(img{1},'sint32');
hl=dipshow(101,ROIs);
h=dipshow(img);
dipmapping(hl,'labels');
n=1;
points={};
fprintf('Please select the ROIs by consecutive clicks. Finish ROI by double click. Undo p9oint by right-click. Finsih all ROIs by double click on a new ROI.\n')
while (1)
    fprintf('Starting new ROI in class number %d. \n',n)
    try
        [roi, v]=diproi(h,'polygon');
    catch
        if isempty(numClasses) || n == numClasses
            break;
        else
            n=n+1;
            continue;
        end
    end
    ROIs(roi)=n;
    points{n}=v;
    if isempty(numClasses)
        n=n+1;
    end
    hl=dipshow(101,ROIs);
end
if isempty(numClasses)
    n=n-1;
end
fprintf('Thank you for selecting %d ROIs. \n',n)

