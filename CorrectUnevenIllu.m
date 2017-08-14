% corrected=CorrectUnevenIllu(img,tilesize,smaller,BgRange,aKernel,ThreshValue)
% Function for a multiplicative correction of the image by an averaged brightness
% img: Image to correct
% tilesize: Size of a single tile (default: [1024 1024])
% smaller: can be supplied to process only a smaller number of tiles (default: [] which means select the whole image)
% BgRange: A number (for background value) or a range from which to estimate the offset)
% aKernel: convolution kernel to average the correction image
% ThreshValue: The correction factor should be estimated only from foreground. This threshold seperated foreground from background. (default: 30)
%
% Authors:  R. Heintzmann, O. Matrosova/Chernavskaia  (8.2.2013)
% for details see
% see F. B. Legesse, O. Chernavskaia, S. Heuke, T. Bocklitz, T. Meyer,  J. Popp, R. Heintzmann, Seamless Stitching of Tile Scan Microscope Images, Journal of Microscopy 258, 223–232,  DOI:10.1111/jmi.12236, 2015

function corrected=CorrectUnevenIllu(img,tilesize,smaller,BgRange,aKernel,ThreshValue)
if nargin < 2
    tilesize=[1024 1024];
end
if iscell(img)
   tilesize=size(img{1,1});
end
if nargin < 3 || isempty(smaller)
if iscell(img)
    smaller=size(img{1,1}).*size(img)./tilesize;
else
    smaller=size(img)./tilesize;
end
end
if nargin < 4 || isempty(BgRange)
    BgRange={0:size(img,1)-1,0:size(img,2)-1};
end
if nargin < 5
    aKernel=200;
end
if nargin < 6
    ThreshValue=30;
end

% Estimate Background
if ~iscell(img)
    img=img(0:tilesize(1)*smaller(1)-1,0:tilesize(2)*smaller(2)-1);
end

if prod(size(BgRange)) == 1
    bg=BgRange;
    fprintf('Using user supplied background value of %g\n',bg);
else
    bg=min(medif(img{1,1}(BgRange{1},BgRange{2}),5));
    fprintf('Estimating Background\n',bg);
% bg=min(img(BgRange(1),BgRange(2)));
    fprintf('Background estimated to be %g\n',bg);
end
% Correct background and make a foreground mask


if iscell(img)
    tilemask=newimar(size(img));
    at=newimar(size(img));
    for p=1:size(img,1)
    for q=1:size(img,2)
        at{p,q}=img{p,q}-bg;
        tilemask{p,q} = at{p,q} > ThreshValue;
    end
    end
else
    % Convert to tiles as cell array (imar)
    img=img-bg;
    amask=img > ThreshValue;
    at = detile(img .* amask,[size(img,2)./tilesize(2) size(img,1)./tilesize(1)]);
    tilemask = detile(amask,[size(img,2)./tilesize(2) size(img,1)./tilesize(1)]);
end


masksum=imarfun('imsum',tilemask * 1.0);
mysum=imarfun('imsum',at) ./ masksum;
xbrightness = mean(mysum,[],2);
ybrightness = mean(mysum,[],1);

plot(xbrightness);hold on;plot(ybrightness,'g');title('Brigness Curves');xlabel('Tile Coordinate');xlabel('Mean Intensity');legend({'x','y'});

xcorrection=gaussf(xbrightness,aKernel) / (mean(xbrightness+ybrightness)/2) ;
ycorrection=gaussf(ybrightness,aKernel) / (mean(xbrightness+ybrightness)/2) ;

hold off;plot(xcorrection);hold on;plot(ycorrection,'g');title('Correction Functions');xlabel('Tile Coordinate');xlabel('Correction Factor');legend({'x','y'});

corrimg=1/(xcorrection * ycorrection);

if iscell(img)
    corrected=newimar(size(img));
    for p=1:size(img,1)
    for q=1:size(img,2)
        corrected{p,q}=at{p,q} .* corrimg;
        corrected{p,q}(corrected{p,q}<0)=0;
    end
    end
else
    corrected=img .* repmat(corrimg,smaller);
end
