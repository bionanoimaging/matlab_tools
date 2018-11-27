% function polImg=CartImg2Pol(img) : tranfers an image on cartesian coordinates to polar coordinates. The radius always steps in distance of one.

function myrad=CartImg2Pol(img,asize,acenter,rrange)
global cuda_enabled;
cu=exist('cuda_enabled','var') & ~isempty(cuda_enabled);
if cu
    disableCuda();
    if isa(img,'cuda')
        img=dip_image_force(img);
    end
end

maxrad=floor(min(size(img)/2));

if nargin < 2 || isempty(asize)
    asize=[400,maxrad];
end
if nargin < 3 || isempty(acenter)
    acenter=MidPos2D(img);
end
if nargin < 3 || isempty(rrange)
    rrange=[0:1:asize(2)];
end

if norm(acenter(1:2)-MidPos2D(img))==0
    ang=phiphi(img);
    myr=rr(img);
else
    px=xx(img,'corner')-acenter(1);
    py=yy(img,'corner')-acenter(2);
    ang=atan2(py,px);
    myr=sqrt(px.^2+py.^2);
end

ascale=1.0*maxrad/2/pi;

% Grid of Destination coordinates
[xq,yq] = meshgrid(ascale*[-pi:2*pi/asize(1):pi], rrange);
% Map onto Grid of Source coordinates
Vq=griddata(ascale*double(ang(:)),double(myr(:))',double(img(:)),xq,yq,'cubic');

%Vq = TriScatteredInterp(double(ang(:))',double(myr(:))',double(img(:))')
%Vq=griddata([double(ang(:));double(myr(:))],double(img(:)));
myrad=dip_image(Vq);
myrad(isnan(myrad))=0;

ang(ang>0)=ang(ang>0)-2*pi;

% [xq,yq] = meshgrid(ascale*[-pi:2*pi/asize(1):pi], rrange);
% Map onto Grid of Source coordinates
Vq=griddata(ascale*double(ang(:)),double(myr(:))',double(img(:)),xq,yq,'cubic');

%Vq = TriScatteredInterp(double(ang(:))',double(myr(:))',double(img(:))')
%Vq=griddata([double(ang(:));double(myr(:))],double(img(:)));
lefthalf=dip_image(Vq);
myrad(xx(myrad)<-1)=lefthalf(xx(myrad)<-1);
myrad(isnan(myrad))=0;
myrad(:,0)=img(floor(size(img,1)/2),floor(size(img,2)/2));  % fix the line at radius zero

if cu
    enableCuda();
end
