% img=squeezedim(img,tosqueeze) : squeezes out dimensions via reshape
% img : image to squeeze
% tosqueeze : a dimension to get rid of. It needs to be singleton, meaning of size 1.
%
% example:
% a =readim();
% b = reshape(a,[1 256 1 1 256]);
% c = squeezedim(b,[1,3]);  % remove dims 1  and 3
function img=squeezedim(img,tosqueeze)
sz = size(img);
tosqueeze(tosqueeze > length(sz)) = [];  % these are ignored as they are already squeezed out
if any(sz(tosqueeze) ~= 1)
    error('squeezedim needs all singleton sizes to squeeze out.')
end
nsz = sz;
nsz(tosqueeze) = [];
img = reshape(img,nsz);
