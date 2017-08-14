%% Creates a movie from a series of 2D RGB images
% It will merge forward and backward directions.
% 
% Synopsis: out=addReverseMovie(filenameOUT, filenameIN, format)
% filenameOUT: filename for the created .avi file
% filenameIN: Name of the first image to read
% format: Format of the individual frames.
% Same as if calling readtimeseries(filenameIN,format)
%
% out is a dip-image array containing the RGB movie
%
% Aurelie, January 2014

function out=addReverseMovie(filenameOUT, filenameIN, format)

% a=readtimeseries('Animation.00','jpg')
a=readtimeseries(filenameIN,format); %a{1} is 4D, the third dimension contains the R,G and B channels.
% a{1}==a{2}; %They are the same

channel1=squeeze(a{1}(:,:,0,:)); %Red channel
channel1F=cat(3,channel1,flipdim(channel1,3)); %Combine forward and backward

channel2=squeeze(a{1}(:,:,1,:)); %Green channel
channel2F=cat(3,channel2,flipdim(channel2,3)); %Combine forward and backward

channel3=squeeze(a{1}(:,:,2,:)); %Blue channel
channel3F=cat(3,channel3,flipdim(channel3,3)); %Combine forward and backward

out=joinchannels('RGB',channel1F,channel2F,channel3F); %Array of images

tmp=cat(4,out{1},out{2},out{3}); %size(tmp)   320   240   200     3
out2=permute(tmp,[1 2 4 3]); %size(out2)   320   240     3   200
filename=[filenameOUT '.avi'];
if (0)
    writeavi(out,'test'); %Apparently, this is soon not going to be supported anymore
else
    vidObj = VideoWriter(filename); %Constructs a video object
    open(vidObj); %Opens (or creates) video file
    writeVideo(vidObj,uint8(out2)); %writes the sequence of images to the video file
    % Second argument must be one of the following classes: double, single,
    % uint8 (no dipimage). It can also not be an array of images (which is the case of out). 
    close(vidObj); % Closes the file.
end
