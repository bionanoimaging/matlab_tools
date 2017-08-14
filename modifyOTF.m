% res=modifyOTF(otf)  : Will apply the pixel form factor effect to an OTF. Assumes square shaped pixels with 100% fill factor
function res=modifyOTF(otf)  
    res=otf .* sinc(xx(otf,'freq')*pi) .* sinc(yy(otf,'freq')*pi);
    