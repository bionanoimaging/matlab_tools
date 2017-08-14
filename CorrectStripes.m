function res=CorrectStripes(data,asigma,dividingline)
if nargin < 2
    asigma=2;
end
if nargin < 3
    dividingline=floor(size(data,2)/2-1);
end
data1=data(:,0:dividingline,:);
data2=data(:,dividingline+1:end,:);
pdata1=mean(data1,[],[2 3]);
pdata2=mean(data2,[],[2 3]);
aratio1=pdata1/gaussf(pdata1,asigma);
aratio2=pdata2/gaussf(pdata2,asigma);
res1=data1/aratio1;
res2=data2/aratio2;
res=cat(2,res1,res2);

