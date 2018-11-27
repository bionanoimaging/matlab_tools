% res=shearArbitrary(img,shearamount,ShearDir,ShearOrtho) : Shears a dataset along one coordinate, in dependence of another coordinate
% img : input image to shear. Size can be 3D or 4D.
% shearamount: amount per ShearOrtho direction to move along direction ShearDir
%

function res=shearArbitrary(img,shearamount,shearDir,shearOrtho)
% if nargin < 3 
%     doExpand=0;
% end
% if doExpand
%     newsize2d=ceil(size2d(img) + abs(shearamount * size(img,3)));
%     img=extract(img,[newsize2d size(img,3) size(img,4)]);  % will keep other dimensions as they are
% end

% shearamount=0.37;
if nargin < 3
    shearDir=3;
end
if nargin < 4
    shearOrtho=4;
end
FTVec=zeros(1,ndims(img)); FTVec(shearDir)=1;
fi=dip_fouriertransform(img,'forward',FTVec);
kx=shearamount*ramp(size(img),shearDir,'freq');
myshear = exp((1i*2*pi)*kx*ramp(size(img),shearOrtho));
res=real(dip_fouriertransform(fi .* myshear,'inverse',FTVec));

