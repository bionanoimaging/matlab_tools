%CAL_READNOISE   Calculates the read noise/gain of a CCD
%
% SYNOPSIS:
%  out = cal_readnoise(in, bg, N, rg, sh)
%
% PARAMETERS:
%  in: input image (a stack of images -at least 10- containing repeated measurements 
%      of the same object); best if all possible intensities are present
%      in the image (create a intensity gradient, e.g. out-of-focus edge) 
%  bg: background image (better a series of background images)
%  N:  number of intensity bins (default 100) x coord is center of a bin
%  rg: fit range, e.g. [1e3 1e4] (default -1, the whole range of intensities)
%  sh: perform shift alignment before evaluation
%
%  out(1) n(e-) variance in the background
%  out(2) the 1/slope of the fit in [photons per ADU]: conversion ADU = detected #photons * slope
%  out(3) n(e-) readnoise fit; sqrt(of y-axis intersection) / slope
%  out(4) mean background
%
%  The gain is obtained by fitting a line to the variance vs intensity plot
%  either weighted or unweighted dependent on which toolboxes are avaibable
%
% LITERATURE:
%  L.J. van Vliet, D. Sudar and I.T. Young in Cell Biology volume III, J.E. Celis (eds.)
%  Digital Fluorescence Imaging Using Cooled CCD Array Cameras, pp.109-120, 1998

% (C) Copyright 2004           Department of Molecular Biology
%     All rights reserved      Max-Planck-Institute for Biophysical Chemistry
%                              Am Fassberg 11, 37077 G"ottingen
%                              Germany
%
% Bernd Rieger, Rainer Heintzmann August 2004.
% added shift correction (Dec 2004)
% bin(0)=bg and weighted fit (BR, KL Dec 2004)
% added 3D background image, some more output (BR,RH July 2007)
% RH 17.04.2016 corrected error in estimating bin brightness. Now from mean value

function out = cal_readnoise(varargin)

d = struct('menu','Statistics',...
   'display','Read noise calibration of a CCD',...
   'inparams',struct('name',       {'a','abg','N','rg','sh'},...
         'description',{'Input image','Background image',...
            'Number of intensity bins','fit range [ADU]','Shift correction'},...
         'type',       {'image','image','array','array','boolean'},...
         'dim_check',  {0,0,0,[1,2],0},...
         'range_check',{[],[],'N+','R',[]},...
         'required',   {1,1,0,0,0},...
         'default',    {'a','b',100,-1,0}...
        ),...
   'outparams',struct('name',{'out'},...
                      'description',{'Output image'},...
                      'type',{'image'}...
                      )...
  );
if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
      out = d;
      return
   end
end
try
   [a,abg,N,rg,sh] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error(firsterr)
   end
end
if ndims(a)~=3
	error('Input stack should be 3D, series of same measurements.');
end

if sh
   [tmp,s]=correctshift(a);
   m=max(max(abs(s)));
   sz=size(a{1});
   a = cut(tmp,[sz(1:2)-ceil(m) sz(3)]);
   fprintf(' Recropping the image to size %d %d %d\n',size(a))
   clear tmp
end

%-------------background correction------------------------
if ndims(abg)==3
   bg = squeeze(mean(abg,[],3));  % calculate a mean bg
   t = a-repmat(bg,[1 1 size(a,3)]); % subtract bg pixel by pixel from data
else
   bg = mean(abg);
   t = a-bg;
end
clear a 
% no = mean(sum(t(:,:,0),[],[1 2]));  % average mean of the first image
no = mean(t);  % average mean of the first image.  Changed by RH, 17.04.2016. This gives a better estimate

%------------normalization of imput data--------------------
% for ii=1:size(t,3)-1
%    t(:,:,ii) = t(:,:,ii)* no/sum(t(:,:,ii));
% end
t=no * t./mean(t,[],[1 2]);  % mormalize for bleaching and intensity fluctuations. Changed RH 17.04.2016

mz = squeeze(mean(t,[],3));    %mean in ADU
v =  squeeze(var(t,[],3));     %variance in ADU

N = N-1; %this leads to N bins
rmz = dip_image(floor(mz/max(mz)*N+.5)+1,'sint32');%quantize the ADU into N bins, use this as label for measure
% The 2 lines below was introduced by RH 17.04.2016
msr = measure(rmz,mz,{'Mean'},[],1,0,0);  % use the mean in the histogram bin to get the x-axis 
tp1 =msr.Mean;  % estimated from mid point of bin (hopefully)
msr = measure(rmz,v,{'Mean','Size'},[],1,0,0);  % get the mean of the variance for the y-axis. "Size" is to know the precision for the RWLS algorithm
% tp1 =[msr.ID] .*max(mz)/N;  % estimated from mid point of bin (hopefully). Commented out RH 17.04.2016
BgVar=double(var(abg));

if ndims(abg)==3
   BgVar=mean(var(abg,[],3));
   mbg=mean(abg,[],3);
   BgVarMean=var(mbg);
   HotImg=mbg-mean(mbg)>3*sqrt(BgVar);
   HotPixels=sum(HotImg);
end

%x = [0 tp1]; %do not include background image in fit as suggested by Rainer
%y = [BgVar msr.Mean];
x = tp1;
y = msr.Mean;
numy = msr.Size;  % how many variance estimates were averaged together.
bgMean=mean(bg);

h=figure;
plot(x+bgMean,y,'x'); hold on
plot(0+bgMean, BgVar,'o'); 
plot([bgMean bgMean], [0 BgVar],'g');
plot([0 bgMean], [BgVar BgVar],'g');
hold off
xlabel('Intensity [ADU]','FontSize',18)
ylabel('Variance [ADU^2]','FontSize',18)
ax = axis;
axis([0 ax(2) 0 ax(4)])
g=gca;
set(g,'FontSize',12)

if ~isempty(rg)
    rg = rg - bgMean;
end
%------------fit-----------------------------------------------
if exist('fit')
   fprintf(' Weigthed fit\n');
   %get the weights (number of pixel per bin)
   msr2 = measure(rmz,rmz,{'Sum'},[],1,0,0);
   w1 = msr2.Sum./msr2.id;
   %szbg = size(abg);
   sz3 = size(t,3);
   w = [w1.*sz3]; %number of pixel that are in that bin
   w = w./sum(w);

   opts = fitoptions;
   if any(rg<0)
      opts = fitoptions(opts,'Weights',w);
      f = fit(x',y','poly1',opts);   
   else
      mask = x > rg(1) & x <rg(2);
      opts = fitoptions(opts,'Weights',w(mask));
      f = fit(x(mask)',y(mask)','poly1',opts);   
   end
   out(1) = f.p2;
   out(2) = f.p1;   
   hold on
   plot([0 x]+bgMean,[out(1) x*out(2)+out(1)],'r');
   hold off
   
   figure;%plot weigths also
   plot(x+bgMean,y,'kx');hold on
   plot([0 x]+bgMean,[out(1) x*out(2)+out(1)],'r');hold off
   xlabel('Intensity [ADU]','FontSize',18)
   ylabel('Variance [ADU^2]','FontSize',18)
   ax = axis;
   axis([0 ax(2) 0 ax(4)]);
   ax1=gca;
   set(ax1,'FontSize',12)   
   ax2 = axes('Position',get(ax1,'Position'),...
      'XAxisLocation','top','YAxisLocation','right',...
      'Color','none','XColor','b','YColor','b');   
   ylabel('relative weight','FontSize',16,'Parent',ax2)   
   if any(rg<0)
      line(x,w);
   else
      line(x(mask),w(mask));
   end  
   axis(ax2,[0 ax(2) 0 .5])
   set(ax2,'YScale','log');
   
   hold off
elseif exist('RWLSPoisson')
   fprintf(' RWLSPoisson fit (Poisson weigthed)\n');
   if any(rg<0)
      [f1,f2,vv]=RWLSPoisson(x,y,numy);
      out=[f1,f2];
   else
      mask = x > rg(1) & x <rg(2);
      [f1,f2,vv]=RWLSPoisson(x(mask),y(mask),numy(mask));
      out=[f1,f2];
   end
   hold on
   fity=[out(1) x*out(2)+out(1)];
   myStd=sqrt([0 vv]);  % add an empty datapoint and convert variance to StdDev.
   % errorbar([0 x]+bgMean,fity,myStd,myStd,'r');
   plot([0 x]+bgMean,fity,'r');
   if any(rg<0)
   else
       x=x(mask);
       fity=[0 fity(mask)];
   end
   plot([0 x]+bgMean,fity-myStd,'r:');
   plot([0 x]+bgMean,fity+myStd,'r:');
   hold off
elseif exist('robustfit')
   fprintf(' robust fit (unweigthed)\n');
   if any(rg<0)
      [out,stats]=robustfit(x,y);
      out=[f1,f2];
   else
      mask = x > rg(1) & x <rg(2);
      [out,stats]=robustfit(x(mask),y(mask));
      out=[f1,f2];
   end
   hold on
   plot([0 x]+bgMean,[out(1) x*out(2)+out(1)],'r');
   hold off
elseif exist('dip_linefit')
   fprintf('dip_line fit (unweigthed)\n');
   if rg<0
      out=dip_linefit(x,y,15);
   else
      mask = x > rg(1) & x <rg(2);
      out=dip_linefit(x(mask),y(mask),15);
   end
   tmp=out(1);out(1)=out(2);out(2)=tmp;
   hold on
   plot(x+bgMean,x*out(2)+out(1),'r');
   hold off
else
   warning('No fitting functionality found.');
   out = -1;
   return   
end
figure(h)

%out(1) = y(1);%readnoise given by variation in background image 
%out(1) = sqrt(out(1))/out(2);

tmp = out;
out(1) = sqrt(BgVar)/tmp(2); % Readnoise background in electrons
out(2) = 1/tmp(2); % gain, to make it photons per ADU
out(3) = sqrt(tmp(1))/tmp(2); % Readnoise fit in electrons
out(4) = bgMean; %offset

if ndims(abg)==3
   s=sprintf('offset: %0.4g ADU\ngain: %0.2g [e^-/ADU],\nreadnoise (Fit): %0.3g e^- RMS,\nreadnoise (Bg): %0.3g e^- RMS,\nfixed pattern noise (Bg): %0.3g e^- RMS\nHot pixels: %3d',...
      out(4),out(2),out(3),out(1),out(2)*sqrt(var(mean(abg,[],3))),HotPixels);
   ax = axis;
   axis([0 ax(2) 0 ax(4)])
   text(ax(2)/10,ax(4)*0.8,s,'FontSize',11)
   fprintf('\nHot Pixels %d\n',HotPixels);
   fprintf('Fixed pattern noise %f\n',out(2)*sqrt(var(mean(abg,[],3))));
   
else
   s=sprintf('gain: %0.2g [e^-/ADU]\nreadnoise (Fit): %0.3g e^- RMS,\nreadnoise (Bg): %0.2g e^- RMS\n',...
   out(2),out(3),out(1));
   h=axis;
   text(x(4),0.8*h(4),s,'FontSize',14)
end
