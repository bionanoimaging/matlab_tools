% function S=AffineSequenceFromMatrix(M) : computes the parameters for the Fourier-shear sequence based on an arbitrary
% Affine transform matrix
% M : Affine transformation matrix
% S : series of transformation steps S1-S9
% S : coefficient vector: 
% S1-S9 shear and scale components
% S10-S12: shift components
% data : data to transform onto new grid
% The shear components have the following meaning:
% Oder of Matrices:  Sm9 * Sm8 * Sm7 * Sm
% Sm1
% 1 S1 S2 
% 0 1  0  
% 0 0  1  
% Sm2
% 1 0  0 
% 0 1 S3  
% 0 0  1  
% Sm3
% S4 0 0 
% 0 S5 0  
% 0 0 S6
% Sm4
% 1  0 0 
% S7 1 0  
% 0  0 1
% Sm5
% 1  0 0 
% 0  1 0  
% S8 S9 1

function S=AffineSequenceFromMatrix(M)
S=zeros(1,9);
S(4) = M(1,1);
S(1) = M(1,2) / M(1,1);
S(2) = M(1,3) / M(1,1);
S(7) = M(2,1) / M(1,1);
T1 = M(1,2)*M(2,1) / M(1,1);
S(5) = M(2,2) - T1;
S(3) = (M(2,3) - M(2,1) * S(2)) / S(5);
S(9) = (M(1,1)*M(3,2) - M(3,1)*M(1,2)) / (M(1,1)*M(2,2) - M(1,2)*M(2,1));
S(8) = (M(3,1) -S(9)*M(2,1)) / M(1,1);
S(6) = M(3,3) -M(1,3)*(S(8)+M(2,1)*S(9)/M(1,1))-S(9)*(M(2,3)-M(2,1)*S(2));

