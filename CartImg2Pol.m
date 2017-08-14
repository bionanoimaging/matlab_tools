% function polImg=CartImg2Pol(img) : tranfers an image on cartesian coordinates to polar coordinates. The radius always steps in distance of one.

function myrad=CartImg2Pol(img,asize)
global cuda_enabled;
cu=exist('cuda_enabled','var') & ~isempty(cuda_enabled);
if cu
    disableCuda();
    if isa(img,'cuda')
        img=dip_image_force(img);
    end
end

maxrad=floor(min(size(img)/2));

if nargin < 2
    asize=[400,maxrad];
end

ang=phiphi(img);
myr=rr(img);

ascale=1.0*maxrad/2/pi;

[xq,yq] = meshgrid(ascale*[-pi:2*pi/asize(1):pi], [0:1:asize(2)]);
Vq=griddata(ascale*double(ang(:)),double(myr(:))',double(img(:)),xq,yq,'cubic');

%Vq = TriScatteredInterp(double(ang(:))',double(myr(:))',double(img(:))')
%Vq=griddata([double(ang(:));double(myr(:))],double(img(:)));
myrad=dip_image(Vq);
myrad(isnan(myrad))=0;

ang(ang>0)=ang(ang>0)-2*pi;

[xq,yq] = meshgrid(ascale*[-pi:2*pi/asize(1):pi], [0:1:asize(2)]);
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
