function out=myRft(myinput)

% b=ft(readim*extract(a,size(readim))). %what I used to test (a was just a
% grating to make the peaks visible)
b=ft(myinput); % dip_image Fourier Transform: whole plane

midX=floor(size(b,1)/2); %middle
midY=floor(size(b,2)/2); %middle

% take the bottom half b(:,midY:end) (X and Y are inverted in dip image). 
% Put the right (midX:end) bottom quarter on the left, 
% and the left (0:midX) bottom quarter on the right.
out=cat(1,b(midX:end,midY:end), b(0:midX,midY:end));