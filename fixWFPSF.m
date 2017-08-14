function h=fixWFPSF(h,arange)
if nargin<2
    arange=0.35;
end
midz=floor(size(h,3)/2);
dz=floor(size(h,3)*arange/2);
exh=h(:,:,midz-dz:midz+dz);
anorm=mean(exh(:,:,0)+exh(:,:,end))/2;
% h(:,:,midz-dz:midz+dz)=h*mean(exh)./mean(exh,[],[1 2]); %Rainer's
% suggestion but dimensions dismatch error
h(:,:,midz-dz:midz+dz)=exh.*anorm./mean(exh,[],[1 2]); %I think that this is what he meant
