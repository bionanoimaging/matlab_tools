% function res=affineFourier(data,M) : performs an affine transform of a 3D dataset by using Fourier-transforms


function [data,S]=affineFourier(data,M, myshift) 
if nargin < 3
    myshift=[0 0 0];
end
if nargin < 2
    al=-pi*45/180;
    M=[cos(al) sin(al) 0;
   -sin(al)  cos(al) 0;  
    0      0  1 ;];
end

allpermutations={[1 2 3],[1 3 2],[2 3 1],[2 1 3],[3 2 1],[3 1 2]};
% S=AffineSequenceFromMatrix(M);  % Transform the affine transform matrix M into the sequence of shear and scale operations
bvval=Inf;
bper=[];bind=0;
bM=[];
for p=1:length(allpermutations)
    myPerm = eye(3);
    myPerm = myPerm(:,allpermutations{p});
    Mper=M * inv(myPerm);
    Sper=AffineSequenceFromMatrix(Mper);  % Transform the affine transform matrix M into the sequence of shear and scale operations
    vval=var(abs(Sper(4:6)));
    if vval < bvval
        bvval = vval;
        bper=allpermutations{p};
        bind=p;
        S=Sper;
        bM=Mper;
        bPerm=myPerm;
    end
%    fprintf('permutationh %d,%d,%d var %g, bestvar %g\n',allpermutations{p},vval,bvval);
end
% calculate all permutations
if bind > 1
    data=permute(data,bper);
end
S=[S myshift];

%%
% sxy=S(1);sxz=S(2);syz=S(3);szx=S(7);szy=S(8);szx2=S(9);

mykxx= xx([size(data,1) 1 1],'freq');mykyy=yy([1 size(data,2) 1],'freq');mykzz=zz([1 1 size(data,3)],'freq');
myxx2= xx([size(data,1) 1 1]);myyy2=yy([1 size(data,2) 1]);myzz2=zz([1 1 size(data,3)]);

rx=exp(-1i*2*pi*mykxx * (S(2) * myzz2 + S(1) * myyy2)); % beam-shear along x in dependence of y an z
rxy=exp(-1i*2*pi*mykyy * (S(3) * myzz2)); % beam shear along y in dependence of z

sxyz=exp(-1i*2*pi*(mykxx * S(10) + mykyy * S(11) +mykzz * S(12))); % shift along x,y and z

% shearings always needs untransformed directions
data=dip_fouriertransform(data,'forward',[1 0 0])*rx;  % x transformed, yz untransformed
data=dip_fouriertransform(data,'forward',[0 1 0])*rxy;  % xy transformed, z untransformed
data=dip_fouriertransform(data,'forward',[0 0 1]);     % xyz transformed

newsize=round(abs(size(data) .* S(4:6)));
data=FourierFlipDim(data,S(4:6) < 0);

% newsize=floor(abs(size(data) .* S(4:6) / 2))*2+1;
% % here scalings can be done
% 
% for d=1:3
%     if (S(3+d) < 0)
%         data=flipdim(data,d);  % accounts for negative scalings, which are not covered by extract
%     end    
% end

data = extract(data .* sxyz,newsize);  % shift and rescale

mykxx= xx([newsize(1) 1 1],'freq');mykyy=yy([1 newsize(2) 1],'freq');mykzz=zz([1 1 newsize(3)],'freq');
myxx2= xx([newsize(1) 1 1]);myyy2=yy([1 newsize(2) 1]);myzz2=zz([1 1 newsize(3)]);
ryz=exp(-1i*2*pi*mykyy * (S(7) * myxx2));
rz=exp(-1i*2*pi*mykzz * (S(9) * myyy2 + S(8) * myxx2));
% now transform backward in a different order
data=dip_fouriertransform(data,'inverse',[1 0 0])*ryz;  % yz transformed, x untransformed
data=dip_fouriertransform(data,'inverse',[0 1 0])*rz;   % z transformed, xy untransformed
data=dip_fouriertransform(data,'inverse',[0 0 1]);     

 fprintf('imaginary ratio is %g\n',max(abs(imag(data))) ./ max(abs(real(data))));
data=real(data);
%if bind > 1
%    data=permute(data,bper);
%end

%% check for the correctness of the calculation
if (0)
    A5 = [1 0 0;0 1 0; S(8) S(9) 1];
    A4 = [1 0 0;S(7) 1 0;0 0 1];
    A3 = [S(4) 0 0;0 S(5) 0; 0 0 S(6)];
    A2 = [1 0 0;0 1 S(3);0 0 1];
    A1 = [1 S(1) S(2);0 1 0; 0 0 1];    
    MyM = A5 * A4 * A3 * A2 * A1;    % should equal M * inv(bPerm) = bM
end

