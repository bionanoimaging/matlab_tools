% h=plotWithAxis(data,stepsize,zeropos,mylabelX,mylabelY,otherargs): easy to use plot function with axis labels and axis ticks
% stepsize : X stepsize per data point e.g. in nm
% zeropos : index to the data where the zero coordinate should be placed
% mylabelX : label for the xaxis (default: 'Position [nm]'
% mylabelY : label for the yaxis (default: 'Intensity [a.u.]'
% otherargs : arguments to pass to the plot function. Have to be enclosed in {}
%
% Example
% plotWithAxis(sin([1:0.1:10]),0.1,5,'Distance [nm]','Intensity [a.u.]')
% saveas(gcf,'TestFigure.png')


function h=plotWithAxis(data,stepsize,zeropos,mylabelX,mylabelY,otherargs)
if nargin < 6
    otherargs={};
end
if nargin < 5
    mylabelY='Intensity [a.u.]';
    if nargin == 4 && iscell(mylabelX)
        otherargs=mylabelX;
        mylabelX='Position [nm]';
    end
end
if nargin < 4
    mylabelX='Position [nm]';
    if nargin == 3 && iscell(zeropos)
        otherargs=zeropos;
        zeropos=0;
    end
end
if nargin < 3
    zeropos=0;
end


myxaxis = stepsize*(0:prod(size(data))-1)-zeropos;

h=plot(myxaxis,data,otherargs{:});
% set(h,'Fontsize',16);
xlabel(mylabelX,'Fontsize',16);
ylabel(mylabelY,'Fontsize',16);
