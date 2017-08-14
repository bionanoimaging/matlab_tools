% GaussArgPermuteXY(myvec,numdims) : Permutes (exchanges) the X and Y coordinates in the arguments to the MultiGaussMSEDeriv function
% arguments are: myvec=[offset, sx,sy, ..., px1,py1, ..., wx1,wy1, ...,int1, px2,py2, .., wx2,wy2, ..., int2, ...]; 
function resvec=GaussArgPermuteXY(myvec,numdims)
if numdims<2
    resvec=myvec;
    return
end
skipdims=(2*numdims+1);
offset=1+numdims;  % backgournd + slopes
numSpots=(length(myvec)-offset)/skipdims;
if ~(numSpots == floor(numSpots))
    error('GaussArgPermuteXY: input vector has wrong length');
end
resvec=myvec;
resvec(offset+1:skipdims:offset+1+skipdims*(numSpots-1))=myvec(offset+2:skipdims:offset+2+numdims*(numSpots-1)); % px and py
resvec(offset+2:skipdims:offset+2+skipdims*(numSpots-1))=myvec(offset+1:skipdims:offset+1+numdims*(numSpots-1));

resvec(offset+1+numdims:skipdims:offset+1+numdims+skipdims*(numSpots-1))=myvec(offset+2+numdims:skipdims:offset+2+numdims+numdims*(numSpots-1)); % wx and wy
resvec(offset+2+numdims:skipdims:offset+2+numdims+skipdims*(numSpots-1))=myvec(offset+1+numdims:skipdims:offset+1+numdims+numdims*(numSpots-1));

resvec(2)=myvec(3); % sx and sy
resvec(3)=myvec(2);
