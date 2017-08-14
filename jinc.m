% res=jinc(mysize,myscales)  : Caculates a bessel(1,2*pi*radius)/radius   = jinc function, which describes the Airy disc in scalar low NA approximation
% example:
% pixelSize=203;  % half the Nyquist freq
% lambda=488/pixelSize;
% na=0.3;
% AbbeLimit=lambda/na;   % coherent Abbe limit, central illumination, not incoherent
% ftradius=1/AbbeLimit*mysize;
% myscales=ftradius./mysize; % [100 100]/(488/0.3);
% res=jinc([256 256],myscales);  % Airy disc

function res=jinc(mysize,myscales)
if nargin < 1
    mysize=[256 256];
end
if nargin < 2
    pixelSize=203;  % half the Nyquist freq
    lambda=488/pixelSize;
    na=0.3;
    AbbeLimit=lambda/na;   % coherent Abbe limit, central illumination, not incoherent
    ftradius=1/AbbeLimit*mysize;
    myscales=ftradius./mysize; % [100 100]/(488/0.3);
end
myradius=pi*rrscale(mysize,myscales);
res=besselj(1,2*myradius) / (myradius);
res(MidPosX(res),MidPosY(res))=1;
