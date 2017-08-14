% [sumBds,mybds,sumBds2,mybds2]=MeanFromCoord(myres,beadsAt,s,myres2) : Extract multiple ROIs at given coordinates, align them and calculate the sum
function [sumBds,mybds,sumBds2,mybds2]=MeanFromCoord(myres,beadsAt,s,myres2)

if numel(s) < size(myres)
    s(size(s):size(s))=s(1);
    fprintf('Warning: Not all sizes for the extraction ROI were given. Using the first size for remaining dimensions');
end

mybds=newim([s size(beadsAt,1)]);
mybds2=newim([s size(beadsAt,1)]);

for b=1:size(beadsAt,1)
    mypos=beadsAt(b,:);
    mybead=extract(squeeze(myres),s,mypos);
    mybds(:,:,b-1)=DampEdge(mybead,0.4,2,1,2);
    if nargin > 3
        mybead2=extract(squeeze(myres2),s,mypos);
        mybds2(:,:,b-1)=DampEdge(mybead2,0.4,2,1,2);
    end
end

ref=squeeze(mybds(:,:,0));  % The first click is the reference
for b=1:size(beadsAt,1)
    si=findshift(squeeze(mybds(:,:,b-1)),ref,'iter');
    mybds(:,:,b-1)=shift(squeeze(mybds(:,:,b-1)),-si);
    if nargin > 3
        mybds2(:,:,b-1)=shift(squeeze(mybds2(:,:,b-1)),-si);
    end
end

sumBds=squeeze(sum(mybds,[],3));
if nargin > 3
    sumBds2=squeeze(sum(mybds2,[],3));
else
    sumBds2=sumBds;
    mybds2=mybds;
end
