% in input image
% rg display range [0 0.1]
% cr where to pu the legend

% Bernd Rieger, 2004.
% changed behaviour of lh, rh too have the low value on the left, Jan 2005

function [out,cr,sz2]=place_lut(in,rg,cr)
if nargin == 1 & ischar(in) & strcmp(in,'DIP_GetParamList') % Avoid being in menu
      out = struct('menu','none');
      return
end

if nargin == 2
   cr = 'auto';
end
in(in<rg(1))=rg(1);
in(in>rg(2))=rg(2);
sz = size(in);

if strcmp(cr,'auto')
   bg = in;
   bg(bg>0)=1; %background by definition of range
   msr = measure(~bg,[],{'Size','Center'},[],2,10,0);
   [i,ii]=max(msr.Size);
   c=msr(ii).Center;
   sz2=sz/2;
   if     sz2(1) > c(1) & sz(2) > c(2) 
      cr = 'lh';
   elseif sz2(1) < c(1) & sz(2) > c(2)
      cr = 'rh';
   elseif sz2(1) < c(1) & sz(2) < c(2)
      cr = 'lv';
   elseif sz2(1) > c(1) & sz(2) < c(2)
      cr = 'rv';
   end   
end

switch cr
   %left right horizontal vertical    
   case 'lh' %oben links
      %tmp = fliplr(xx(round([sz(1)/3 sz(2)/20])));
      tmp = (xx(round([sz(1)/3 sz(2)/20])));
      sz2 = size(tmp);
      cr = [10 10];
   case 'rh' %oben rechts
      tmp = (xx(round([sz(1)/3 sz(2)/20])));
      sz2 = size(tmp);
      cr = [sz(1)-10-sz2(1) 10];
   case 'lv' %unten links
      tmp = flipud(yy(round([sz(1)/20  sz(2)/3])));
      sz2 = size(tmp);
      cr = [10 sz(2)-10-sz2(2)];
   case 'rv' %unten rechts
      tmp = flipud(yy(round([sz(1)/20  sz(2)/3])));
      sz2 = size(tmp);
      cr = [sz(1)-10-sz2(1) sz(2)-10-sz2(2)];   
   otherwise
      error('Unkown placing option');
end
tmp = stretch(tmp,0,100,rg(1),rg(2));
out = in;
out(cr(1):cr(1)+sz2(1)-1,cr(2):cr(2)+sz2(2)-1) = tmp;


