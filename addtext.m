function res=addtext(img,mytxt,fs,brightness)
if nargin < 3 || isempty(fs)
    fs=10;
end
if nargin < 4
    brightness=255;
end
h=figure();
axis tight;axis off;axis ij; axis([0 1 0 1]);
set (h,'Color','white');
t=text(fs/size(img,1)/20,fs/(fs+16)/15,mytxt,'fontsize',fs);
q=getframe();
q=255-dip_image(q.cdata);
q=q/255*brightness;
res=cat(2,img,squeeze(q(0:size(img,1)-1,0:fs+16,0)));
close(h);
