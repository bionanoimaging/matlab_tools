%% Shows the peaks in Fourier space in SIM data based on their auto-correlation
% 
% showFourierPeaks(data,numdirs,mypsf,displayexp)
%
% data: SIM raw data
% numdirs: number of direction
% mypsf: point spread function. If the user enters a number, a theoretical PSF of radius mypsf will be calculated
% displayexp: exponent to tune for better display. Default 1

function showFourierPeaks(data,numdirs,mypsf,displayexp)
if nargin < 4
    displayexp=1;
end
if isnumeric(mypsf)
    h=rr([size(data,1) size(data,2)],'freq')< mypsf;
    h=ft(ift(h).*conj(ift(h)));
else
    h=ft(mypsf);
end

mypsf=real(ift(h));

data=DampEdge(data,0.1,2,1,2);

af=real(dip_fouriertransform(dip_fouriertransform(data,'forward',[1 1 0]) .* h,'inverse',[1 1 0]));

normalizer=ft(mypsf.*mypsf).^displayexp;
% 
% q=dip_fouriertransform(af.*af,'forward',[1 1 0]) ./ (normalizer + 1e-4*max(abs(normalizer))); % calculate the autocorrelation
% 
% qq=abs(q);
% 
% phases=floor(size(data,3)/numdirs);
% for n=0:numdirs-1
%     sum(qq(:,:,n*phases:(n+1)*phases-1),[],3) .^0.3
% end
% 

phases=floor(size(data,3)/numdirs);
for n=0:numdirs-1
    s=mean(af(:,:,n*phases:(n+1)*phases-1),[],3);
    
    q=dip_fouriertransform(s.*(af(:,:,n*phases:(n+1)*phases-1)-s),'forward',[1 1 0]) ./ (normalizer + 1e-4*max(abs(normalizer))); % calculate the autocorrelation
    
    qq=abs(q);
    
    sum(qq,[],3) .^0.3
end


% for n=0:numdirs-1   % Correlate with the extracted first order
%     s=ft(af(:,:,n*phases:(n+1)*phases-1));
%     s=ift(s(:,:,floor(size(s,3)/2)-2));
%     
%     q=dip_fouriertransform(s.*af,'forward',[1 1 0]) ./ (normalizer + 1e-4*max(abs(normalizer))); % calculate the autocorrelation
%     
%     qq=abs(q);
%     
%     phases=floor(size(data,3)/numdirs);
%     sum(qq,[],3) .^0.3
% end
