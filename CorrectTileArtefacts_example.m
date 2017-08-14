%  Example about how to use RemoveTileArtefacts.m and CorrectUnevenIlluminaiton.m
% 
% Authors:  R. Heintzmann, O. Matrosova/Chernavskaia  (8.2.2013)
% for details see
% see F. B. Legesse, O. Chernavskaia, S. Heuke, T. Bocklitz, T. Meyer,  J. Popp, R. Heintzmann, Seamless Stitching of Tile Scan Microscope Images, Journal of Microscopy 258, 223–232,  DOI:10.1111/jmi.12236, 2015

clear all
% a1=readim('KH50CARS128.tif');
cd \\mars\user\rheintzmann\MATLAB\TobiasMeyer\
a1=readim('Data/KH50CARS.tif');  % a 7680x7680 example dataset
% [scale, offset] = pcfo(a1(0:4:end,0:4:end));  % Not needed! Just used here to estimate the offset (where is the zero intensity?)
offset = 0;  % zero intensity level

tilesize=[512 512];
%tilesize=[128 128];
smaller=size(a1) ./ tilesize;
% Determine and Correct the mean uneven illumination
%corrected=CorrectUnevenIllu(a1,tilesize,smaller,offset,3,10);
corrected=CorrectUnevenIllu(a1,tilesize,smaller,offset,8,0);

% tiffwrite('Corrected.tif',corrected);
% writeim(corrected,'Corrected.tif');
a1r=resample(a1,[0.2 0.2])  % just to speed the calculation up
c1r=resample(corrected,[0.2 0.2])

corrected(corrected<1)=1;
%res=RemoveTileArtefacts(corrected,tilesize,0,[],4,1);
res=RemoveTileArtefacts(corrected,tilesize,0,[],50,1);  % correct the full resolution image.

% dipshow(res)  % shows the full result

%res=RemoveTileArtefacts(corrected,tilesize,0,[],50,10);
res_res=resample(res,[0.2 0.2])  % for display purposes only

cat(3,a1r,c1r,res_res)  % just to compare the images visually. Press "n" and "p" to switch between images

%%
%tiffwrite('Result_Big2.tif',res,'yes');
%tiffwrite('Result_Big2resampled.tif',res_res,'yes');

res=RemoveTileArtefacts(corrected,tilesize,0,[],20,1,[3 3],[3 3]);  % Avoids the dark stripes
dipshow(res)  % display the result
%tiffwrite('Result_3.tif',res);
%tiffwrite('Result_3resampled.tif',res_res);
