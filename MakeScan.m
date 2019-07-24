% scan=MakeScan(sz,Sx,Sy) : creats a scanning deltafunction moving along X and Y
% sz : size of image to generate
% Sx : steps to scan in each scanstep along X
% Sy : steps to scan in each scanstep along Y
%
% Example:
% MakeScan([256 256], 20,20)
function scan=MakeScan(sz,Sx,Sy)

N=prod(floor((sz-1)./[Sx,Sy]+1));
scan=newim([sz N]);

xind = 0:Sx:sz(1)-1;  % fast way based on one-D indexing
yind = 0:Sy:sz(2)-1;
allind = xind + (yind*sz(1))'; % list of indices in xy
allind = allind(:) + sz(1)*sz(2)*[0:N-1]';
scan(allind)=1;

% n=0;
% for y=0:Sy:sz(2)-1   % old slow way
%     for x=0:Sx:sz(1)-1
%         scan(x,y,n)=1.0;
%         n=n+1;
%     end
%     fprintf('Line %d/%d\n',y,sz(2));
% end
