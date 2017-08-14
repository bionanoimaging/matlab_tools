% res=RemoveTileArtefacts(img,tilesize,adir,smaller,sigma,eps,overlap) : Removed the edges between tiles in a mosaic image
% The algorithm is based on Comparing the neighbouring edges and forcing brightnesses at the corners to the average between the mean edges
% img: Input image
% tilesize: size of a single tile (default: [1024 1024])
% adir: direction to apply correction to (1 : x, 2: y, 0: both) default: 0
% smaller: use only a smaller number of tiles (default: [] = all tiles)
% sigma: smoothing kernelsize for measuring edge brightness. Default: 30
% eps: Regularisation parameter, if the average brightness reaches zero. Default: 1
% overlap: 2D vector of number of pixels that overlap between neighboring tiles
% stitchwidth: 2D vector of number of pixels to account in the overlap region
% skew: 2D vector of number of pixels to account in skewing the tiles (for each x tile it is shifted a skew along y and viceversa)
%
% Authors:  R. Heintzmann, O. Matrosova/Chernavskaia  (8.2.2013)
% for details see
% see F. B. Legesse, O. Chernavskaia, S. Heuke, T. Bocklitz, T. Meyer,  J. Popp, R. Heintzmann, Seamless Stitching of Tile Scan Microscope Images, Journal of Microscopy 258, 223–232,  DOI:10.1111/jmi.12236, 2015



function res=RemoveTileArtefacts(img,tilesize,adir,smaller,sigma,eps,overlap,stitchwidth,skew)
if nargin < 9
    skew=[0 0];
end
if nargin < 8
    stitchwidth=[0 0];
end
if nargin < 7
    overlap=[0 0];
end
if nargin < 2 
   tilesize=[1024 1024];
end
if iscell(img)  || isa(img,'dip_image_array')
   tilesize=size(img{1,1});
end
if nargin < 3
    adir=0;
end
if nargin < 4 || isempty(smaller)
if iscell(img) ||  isa(img,'dip_image_array')
    smaller=size(img{1,1}).*size(img)./tilesize;
else
    smaller=size(img)./tilesize;
end
end
if nargin < 5
    sigma=200;
end
if nargin < 6
    eps=1;
end
%smaller=[15 10];
%a1=a1(0:tilesize(1)*smaller(1)-1,0:tilesize(1)*smaller(2)-1);

if iscell(img) || isa(img,'dip_image_array')
    at=img;
    if norm(skew) > 0
        fprintf('Korrecting skew ...');
        for nx=1:size(at,1)
            for ny=1:size(at,2)
                at{nx,ny}=circshift(at{nx,ny},[skew(1)*(ny-1),skew(2)*(nx-1)]);
            end
        end
        fprintf('... done\n');
    end
    at=reshape(cat(3,at{:}),[tilesize(1) tilesize(2) smaller(1) smaller(2)]); % makes it 4D
else
    at = detile(img,[size(img,2)./tilesize(2) size(img,1)./tilesize(1)]);  % makes a cell array
    at=reshape(cat(3,at'),[tilesize(1) tilesize(2) smaller(1) smaller(2)]); % makes it 4D
end
at=reshape(cat(3,at'),[tilesize(1) tilesize(2) smaller(1) smaller(2)]); % makes it 4D

%% Now this is a 4D dataset, dimension 3: to the right, dimension 4 to the bottom


% all_top=gaussf(at(:,0,:,:),[sigma 0 0 0]);
% all_bot=gaussf(at(:,end,:,:),[sigma 0 0 0]);
% all_left=gaussf(at(0,:,:,:),[0 sigma 0 0]);
% all_right=gaussf(at(end,:,:,:),[0 sigma 0 0]);
overlap=floor(overlap/2);
stitchwidth=floor(stitchwidth/2);

% Reduces the area also to what is necessary
all_t=at(overlap(1)-stitchwidth(1):end-overlap(1)+stitchwidth(1),overlap(2)-stitchwidth(2):overlap(2)+stitchwidth(2),:,:);
all_b=at(overlap(1)-stitchwidth(1):end-overlap(1)+stitchwidth(1),end-overlap(2)-stitchwidth(2):end-overlap(2)+stitchwidth(2),:,:);
all_l=at(overlap(1)-stitchwidth(1):overlap(1)+stitchwidth(1),overlap(2)-stitchwidth(2):end-overlap(2)+stitchwidth(2),:,:);
all_r=at(end-overlap(1)-stitchwidth(1):end-overlap(1)+stitchwidth(1),overlap(2)-stitchwidth(2):end-overlap(2)+stitchwidth(2),:,:);

all_top=gaussf(mean(all_t,[],2),[sigma 0 0 0]);
all_bot=gaussf(mean(all_b,[],2),[sigma 0 0 0]);
all_left=gaussf(mean(all_l,[],1),[0 sigma 0 0]);
all_right=gaussf(mean(all_r,[],1),[0 sigma 0 0]);

% determine mean curves on edges
avg_x=cat(4,all_top(:,:,:,0),(all_bot(:,:,:,0:end-1) + all_top(:,:,:,1:end))/2,all_bot(:,:,:,end));
avg_y=cat(3,all_left(:,:,0,:),(all_right(:,:,0:end-1,:) + all_left(:,:,1:end,:))/2,all_right(:,:,end,:));

% determine x and y connection points
connection_x = cat(3,avg_x(overlap(1),:,0,:),(avg_x(end-overlap(1),:,0:end-1,:) + avg_x(overlap(1),:,1:end,:))/2,avg_x(end-overlap(1),:,end,:));
connection_y = cat(4,avg_y(:,overlap(2),:,0),(avg_y(:,end-overlap(2),:,0:end-1) + avg_y(:,overlap(2),:,1:end))/2,avg_y(:,end-overlap(2),:,end));

% corners
corners = (connection_x + connection_y)/2;

%% Now compute correction curves
% stretch endings to corners
x_slope=(xx(avg_x,'corner')-stitchwidth(1))/(size(avg_x,1)-1-2*stitchwidth(1));  % runs from -1 to 1
delta_left=avg_x(stitchwidth(1),:,:,:)-corners(:,:,0:end-1,:);
delta_right=avg_x(end-stitchwidth(1),:,:,:)-corners(:,:,1:end,:);
avg_x = avg_x - ((1-x_slope)*delta_left + x_slope*delta_right);

y_slope=(yy(avg_y,'corner')-stitchwidth(2))/(size(avg_y,2)-1-2*stitchwidth(2));  % runs from -1 to 1
delta_top=avg_y(:,stitchwidth(2),:,:)-corners(:,:,:,0:end-1);
delta_bot=avg_y(:,end-stitchwidth(2),:,:)-corners(:,:,:,1:end);
avg_y = avg_y - ((1-y_slope)*delta_top + y_slope*delta_bot);

% Now compute a correction map
cor_top=abs(avg_x(:,:,:,0:end-1)) ./ (abs(all_top) + eps);
cor_bot=abs(avg_x(:,:,:,1:end)) ./ (abs(all_bot) + eps);
cor_left=abs(avg_y(:,:,0:end-1,:)) ./ (abs(all_left) + eps);
cor_right=abs(avg_y(:,:,1:end,:)) ./ (abs(all_right) + eps);

x2_slope=x_slope(:,:,:,0:end-1);
y2_slope=y_slope(:,:,0:end-1,:);

x2_slope(x2_slope<0)=0;x2_slope(x2_slope>1)=1;
y2_slope(y2_slope<0)=0;y2_slope(y2_slope>1)=1;

% Test it using only the map_x

%%
%resample(corrected,[0.2 0.2])
if adir == 1
    map_x = ((1-x2_slope)*cor_left + x2_slope*cor_right);
    res=reshape(permute(at .* map_x,[1 3 2 4]),size(img));
elseif adir == 2
    map_y = ((1-y2_slope)*cor_top + y2_slope*cor_bot);
    res=reshape(permute(at .* map_y,[1 3 2 4]),size(img));
else
    % p=6;
    p=6;
    % p=12;
    eps=0.5;
    %eps=0.01;
    weight_l = 1/(x2_slope+eps)^p - 1/(1+eps)^p;  % Weighting functions (Sheppards method)
    weight_r = 1/(1-x2_slope+eps)^p  - 1/(1+eps)^p;
    weight_t = 1/(y2_slope+eps)^p  - 1/(1+eps)^p;
    weight_b = 1/(1-y2_slope+eps)^p  - 1/(1+eps)^p;
    fprintf('assembling map ... this may take a while ... ')
    map = (weight_l*cor_left + weight_r*cor_right + weight_t*cor_top + weight_b*cor_bot) ./ (weight_l+weight_r+weight_t+weight_b);
    fprintf('done \n')
    % reshape(permute(map,[1 3 2 4]),size(img))
%% Interpolate the seams
    at=at(overlap(1)-stitchwidth(1):end-overlap(1)+stitchwidth(1),overlap(2)-stitchwidth(2):end-overlap(2)+stitchwidth(2),:,:);
    at=at .* map;
if (1)  % write the weighted averaged data into the tile edge regions
    all_t=at(:,0:2*stitchwidth(2),:,:);;
    all_b=at(:,end-2*stitchwidth(2):end,:,:);
    all_l=at(0:2*stitchwidth(1),:,:,:);
    all_r=at(end-2*stitchwidth(1):end,:,:,:);

    weight_hor=xx([size(all_l,1) size(all_l,2)],'corner')/(size(all_l,1)-1); %linear weight
    mean_hor=all_r(:,:,0:end-1,:).*(1-weight_hor)+all_l(:,:,1:end,:).*weight_hor;
    weight_vert=yy([size(all_t,1) size(all_t,2)],'corner')/(size(all_t,2)-1); %linear weight
    mean_vert=all_b(:,:,:,0:end-1).*(1-weight_vert)+all_t(:,:,:,1:end).*weight_vert;
    
    at(:,0:2*stitchwidth(2),:,1:end)=mean_vert;
    at(:,end-2*stitchwidth(2):end,:,0:end-1)=mean_vert;
    at(0:2*stitchwidth(1),:,1:end,:)=mean_hor;
    at(end-2*stitchwidth(1):end,:,0:end-1,:)=mean_hor;
end
    
    at=at(stitchwidth(1):end-stitchwidth(1),stitchwidth(2):end-stitchwidth(2),:,:);
    res=reshape(permute(at,[1 3 2 4]),[size(at,1)*size(at,3) size(at,2)*size(at,4)]);
end
