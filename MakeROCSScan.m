% illuamp=MakeROCSScan(sz,PupilRadius,NumIllu) : creats a series of azimutally rotating fields, according to Alexander Rohrbach's ROCS technique
% PupilRadius : distance to center the Fouier space of the plane waves
% NumIllu : number of illuminations
% Example:
% MakeROCSScan([256 256], 20,32)
function illuamp=MakeROCSScan(sz,PupilRadius,NumIllu)

krel=PupilRadius./(sz/2);
kang = (zz([1 1 NumIllu],'freq')+0.5)*2*pi;
kvecX=krel(1).*cos(kang);
kvecY=krel(2).*sin(kang);

illuamp=exp(1i*pi*(xx([sz(1) 1 1])*kvecX+yy([1 sz(2) 1])*kvecY));
