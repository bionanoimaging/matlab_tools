% [sumBds,mybds,sumBds2,mybds2]=MeanFromCoord(myres,beadsAt,s,myres2) : Extract multiple ROIs at given coordinates, align them and calculate the sum
function [sumBds,mybds,sumBds2,mybds2,si]=MeanFromCoord(myres,beadsAt,s,myres2,si)
if nargin < 5
    si={};
end
myres=squeeze(myres);  % get rid of an empty 3rd dimension

if numel(s) < ndims(myres)
    s(size(s,2)+1:ndims(myres))=s(1);
    fprintf('Warning: Not all sizes for the extraction ROI were given. Using the first size for remaining dimensions\n');
end

%mybds=newim([s size(beadsAt,1)]);
%mybds2=newim([s size(beadsAt,1)]);
mybds={};
mybds2={};

for b=1:size(beadsAt,1)
    mypos=beadsAt(b,:);
    mybead=extract(squeeze(myres),s,mypos);
    mybds{b}=DampEdge(mybead,0.4,2,1,2);
    if nargin > 3 && ~isempty(myres2)
        mybead2=extract(squeeze(myres2),s,mypos);
        mybds2{b}=DampEdge(mybead2,0.4,2,1,2);
    end
end

ref=squeeze(mybds{1});  % The first click is the reference
for b=1:size(beadsAt,1)
    if numel(si) < b
        si{b}=findshift(squeeze(mybds{b}),ref,'iter');
    end
    mybds{b}=shift(squeeze(mybds{b}),-si{b});
    if nargin > 3  && ~isempty(myres2)
        mybds2{b}=shift(squeeze(mybds2{b}),-si{b});
    end
end

sumBds=squeeze(sum(cat(4,mybds{:}),[],4));
if nargin > 3  && ~isempty(myres2)
    sumBds2=squeeze(sum(cat(4,mybds2{:}),[],4));
else
    sumBds2=sumBds;
    mybds2=mybds;
end
