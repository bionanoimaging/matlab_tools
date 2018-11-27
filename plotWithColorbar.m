%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plotWithColorbar(data, varargin) plots a dip_image input data as follows: 
%
%
% Example:
% plotWithColorbar(readim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotWithColorbar(data, varargin)
%Xaxis, Yaxis, scaleAxis, cutout, center, mylabelX, mylabelY, mytitle, climit
%%
options = struct('Xaxis',(0:size(data,1)-1),'Yaxis',(0:size(data,2)-1),'scaleAxis',[1 1],'cutout',[size(data,1) size(data, 2)],...
    'center',floor(size(data)/2),'mylabelX','X label (a.u.)','mylabelY','Y label (a.u.)','mytitle','Figure Title','climit',[min(min(data)),max(max(data))]);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('GenObj needs paramName/paramValue pairs')
end

for pair = reshape(varargin,2,[]) % pair is {paramName;paramValue}
   inpName = pair{1}; 

   if any(strcmp(inpName,optionNames)) % looks for matching parameter names
      % overwrite parameters
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

for i=1:length(optionNames)
    str=[optionNames{i} '=options.' optionNames{i} ';']; % for example obj=options.obj
    eval(str); % executes obj=options.obj;
end

%% 

%Extract the data of size(cutout) from the given data 
z = double(extract(data, cutout, center));

%Set up plot axis to start from the center of the plot cutout
tmpX = (-floor(cutout(1)/2):floor(cutout(1)/2)-1).*scaleAxis(1);
tmpY = (-floor(cutout(2)/2)+1:floor(cutout(2)/2)).*scaleAxis(2);

%Plot data
surf(tmpX, tmpY, z, 'edgecolor','none', 'facecolor', 'interp');
view(2);
%set(gcf,'renderer','zbuffer');
set(gcf,'renderer','painters'); % to avoide white artefacts in the plot

%Label axes
xlabel(mylabelX,'Fontsize',16);
ylabel(mylabelY,'Fontsize',16);
title(mytitle,'Fontsize',16);

%Add colorbar
caxis(climit)
colorbar
axis tight

