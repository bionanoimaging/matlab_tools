% colFig=showCpxFig(img,myexponent,lut_cr)  : Displayes a complex valued image with brighness and color
% Example:
% showCpxFig(ft(readim),0.3,'lh');
% 

function colFig=showCpxFig(img,myexponent,lut_cr)
if nargin < 2
    myexponent=0.2;
end
mymag=abs(img)^myexponent;
myphase=phase(img);
if nargin >= 3
    [myphase2,cr,sz2] = place_lut(myphase,[-pi pi],lut_cr);
    mymag(cr(1):cr(1)+sz2(1)-1,cr(2):cr(2)+sz2(2)-1)=max(mymag);
%    mymag(myphase2~=myphase)=max(mymag);  % THis is a hack!
    myphase=myphase2;
end
myphase=myphase/2/pi;

R=clipit(myphase);
G=clipit(myphase+1/3);
B=clipit(myphase+2/3);
colFig=joinchannels('RGB',mymag*R,mymag*G,mymag*B);

if nargin >= 3
dipshow(colFig);
fs = 16;%fontsize
txtcol = 'white';

switch lut_cr
   case {'rh','lh'}
     text(cr(1),cr(2)+sz2(2)+fs,'-\pi','Color',txtcol,...
        'FontSize',fs);
     text(cr(1)+sz2(1),cr(2)+sz2(2)+fs,'\pi','Color',txtcol,...
        'FontSize',fs,'HorizontalAlignment','right');
   case 'rv'
      text(cr(1)-sz2(1),cr(2),'\pi','Color',txtcol,...
        'FontSize',fs);
      text(cr(1)-sz2(1),cr(2)+sz2(2),'-\pi','Color',txtcol,...
        'FontSize',fs);
   case 'lv'
      text(cr(1)+sz2(1),cr(2),'\pi','Color',txtcol,...
        'FontSize',fs);
      text(cr(1)+sz2(1),cr(2)+sz2(2),'-\pi','Color',txtcol,...
        'FontSize',fs);
end
end
end

function val=clipit(val)
val=mod(val,1);
val=1-abs(val*3-1.0);
val(val<0)=0;
end