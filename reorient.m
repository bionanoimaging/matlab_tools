function res=reorient(vec,d)
% res=reorient(vec,d) : reorients a 1D dipimage in space
% vec : vector to orient
% d : dimension
    sizes=ones(1,d);
    sizes(d) = prod(size(vec));
    res= reshape(vec,sizes);
    