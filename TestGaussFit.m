a=4+51.2*exp(-(((xx(20,20,10)-3.2)/2.2).^2+((yy(20,20,10)+2.8)/1.7).^2+(zz(20,20,10)/2.6).^2))
b=60.2*exp(-(((xx(20,20,10)-1.2)/2.2).^2+((yy(20,20,10)+0.8)/1.7).^2+((zz(20,20,10)-2)/2.6).^2))

an=noise(a+b,'poisson')
[params,res,fitted,residual]=FitDataNDFast([4 0 0 0 3 -2 0.5 2 2 2 50 1 1 1 2 2 2 50],an,3,300,'mse')
% now lets fix some of the parameters (no slopes)
[params,res,fitted,residual]=FitDataNDFast([4 0 0 0 3 -2 0.5 2 2 2 50 1 1 1 2 2 2 50],an,3,300,'mse',[1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1])

method=1;
for q=1:100
    an=noise(a,'poisson');
    switch method
        case(1)
        [params,res,fitted,residual]=FitDataNDFast([10 0 0 0 3 -3 1 3 3 2 40],an,3,300,'mse',[1 1 1 1 1 1 1 1 1 1 1]);
        case(2)
        [params,res,fitted,residual]=FitDataNDFast([10 0 0 0 3 -3 1 3 3 2 40],an,3,300,'mse',[1 0 0 0 1 1 1 1 1 1 1]);
        case(3)
        [params,res,fitted,residual]=FitDataNDFast([10 0 0 0 3 -3 1 3 3 2 40],an,3,300,'fidiv',[1 0 0 0 1 1 1 1 1 1 1]);
    end
    allparams(q,:)=params;
end
std(allparams(:,5:7))
