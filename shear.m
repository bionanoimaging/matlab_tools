% res=shear(img,shearvec,doExpand) : (de-)shears a dataset by the shear vector (per Z-slice)
% img : input image to shear. Size can be 3D or 4D.
% shearvec : amount of shear in X and Y for each Z plane
% doExpand : boolean, if true, the dataset will be expanded to account for the maximal shear.
%

function res=shear(img,shearvec,doExpand)
if nargin < 3 
    doExpand=0;
end
if doExpand
    newsize2d=ceil(size2d(img) + abs(shearvec * size(img,3)));
    img=extract(img,[newsize2d size(img,3) size(img,4)]);  % will keep other dimensions as they are
end

fi=ft2d(img);
kx=shearvec(1)*xx(size2d(img),'freq');
ky=shearvec(2)*yy(size2d(img),'freq');
myshear = exp((i*2*pi)*(kx+ky)*zz([1 1 size(img,3)]));
res=real(ift2d(fi .* myshear));

