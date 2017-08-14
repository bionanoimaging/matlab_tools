%% Corrects the drift between slices of a dataset
% [res,myshift]=correctZdrift(in,filename,doPlot,alternateAlign)
% in: input dataset (3D)
% filename: Name of the txt file in which the applied shifts are saved. Default: MyShifts.txt (but will overwrite this file if it already exists)
% doPlot: if 1, will plot the shifts with respect to the slice number. Default 1
% alternateAlign: Single plane image to align to. 
% If alternateAlign is NOT given, the alignment will be sequential from plane to next plane.
% This can be very unstable for larger datasets, but has the advantage to also work for Z-stacks with different information in each plane.
% Use a fixed referenc image, if all images look quite similar!
%
% Aurelie, January 2014, Rainer July 2015

function [res,myshift]=correctZdrift(in,filename,doPlot,alternateAlign)

if nargin < 2 || isempty(filename)
    filename=[tempdir 'ResShifts.txt'];
end

if nargin < 3
    doPlot=1;
end

if ~isstr(filename) %filename is not a string
    error('Input filename must be a string');
end

if nargin < 4
    res=kcorrelator(in,{'p',1;'v',filename});
else
    in=cat(3,alternateAlign,in);
    res=kcorrelator(in,{'p',1;'v',filename;'fixplane',0});
    res=res(:,:,1:end);
end
myshift=load(filename);
%midslice=floor(size(myshift,1)/2)+1+mod(size(myshift,1),2);
if nargin < 4
    midslice=ceil(size(myshift,1)/2)+1;
    myshift=[myshift(end:-1:midslice,:);0 0 0;myshift(1:midslice-1,:)];
else
    myshift=-myshift;
end

if doPlot
    plot(myshift(:,1)) %X shift
    hold on %to display on same graph
    plot(myshift(:,2),'r') %Y shift (in red)
    legend('X-shift','Y-shift')
    xlabel('Slice number')
    ylabel('Shift in pixels')
end
 

