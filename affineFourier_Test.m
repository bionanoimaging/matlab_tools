% function affineFourier(M,S,data) : performs an affine transform of a 3D dataset by using Fourier-transforms
% M : affine 3x3 transform matrix
% S : shift vector
% data : data to transform onto new grid

% [U D V]= svd(M);  % D contains the zooms, V contains the 3D rotation

% Now the rotation has to be written as consicutive shearings along X, Y and Z
% A shear in X ist written as
% 1 mxy 0     1 0 mxz  1  0 0
% 0 1  0      0 1 0    mx 1 0
% 0 0  1      0 0 1    0  0 1

%% Make a test object
a=newim([124 124 3]);
ax=abs(xx(a));
ay=abs(yy(a));
az=abs(zz(a));
smin=17; smax=22;
a(ax < smax & ay < smax & az < smax)=3;
a=a* (rro(size2d(a),[15 10]) > 5);
a=gaussf(a,0.7);

%%
% a=readim('Y:\MATLAB\images\resolution_fine.tif')
a=readim('Y:\MATLAB\images\resolution_coarse.tif')
a=cat(3,a,a)
%% try a Fourier-transform with shearings
al=pi*50/180;
rotation3dFourier(a,[al 0 0])
rotation3dFourier(a,[0 al 0])
rotation3dFourier(a,[0 0 al])
rotation3dFourier(a,[al al al])

%% 36 sequential rotations by 10 degrees.
r=a;
rBL=a;
rLA=a;
N=33;
for p=1:N
    p
    deltaAngle=[360/N*pi/180 0 0];
    r=rotation3dFourier(r,deltaAngle);
    rLA=rotation3d(rLA,'direct',deltaAngle,[],'linear');
    rBL=rotation3d(rBL,'direct',deltaAngle,[],'lanczos4');
end

cat(4,a,extract(r,size(a)),rLA,rBL,a)

%%
deltaAngle=[50*pi/180 0 0];
[r1,S1]=rotation3dFourier(a,deltaAngle );
r1LA=rotation3d(a,'direct',deltaAngle,[],'linear');
r1BL=rotation3d(a,'direct',deltaAngle,[],'bspline');
cat(4,extract(r1,size(a)),r1LA,r1BL)

[r2,S2]=rotation3dFourier(r1,-deltaAngle);

cat(4,a,extract(r2,size(a)))

%%
B=50*pi/180;
if (0)
    M=[cos(all) sin(all) 0;
    -sin(all)  cos(all) 0;  
    0      0  1];
else
%  M=[cos(al) 0 sin(al);0      1   0 ; -sin(al)  0 cos(al)];
 M=[1,0,0 ; 0,cos(b),sin(b);0,-sin(b),cos(b)];
end
% M=eye(3);

% S=[0.5 0.4 0.3 1.2 1.2 1.2 0.8 -0.3 0.5  10 10 10];
[res,S]=affineFourier(a,M)

