% GenMask(pointList,radius,asize) : generates a mask based on a list of points
% pointList : List of points
% radius: radius of mask in each point to generate
% asize: size of the mask image to generate
%
% Example:
% GenMask([10,10;120,120;120,17],6,[256 256])

function mymask=GenMask(pointList,radius,asize)
mymask=newim(asize);
mykernel=rr(asize)<radius;
for n=1:size(pointList,1)
    mypos=floor(pointList(n,:)+0.5);
    if numel(mypos) == 2
        if mypos(1) > 0 && mypos(1) < asize(1) && mypos(2) > 0 && mypos(2) < asize(2)
            mymask(mypos(1),mypos(2))=1;
        end
    else
        if mypos(1) > 0 && mypos(1) < asize(1) && mypos(2) > 0 && mypos(2) < asize(2) && mypos(3) > 0 && mypos(3) < asize(3)
            mymask(mypos(1),mypos(2),mypos(3))=1;
        end
    end
end

mymask=convolve(mymask,mykernel) > 0.5;
