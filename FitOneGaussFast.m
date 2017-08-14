% [spots, res, fitted, residual] = FitOneGaussFast(Bg, mystart, Meas) : Fits a Gauss in a non-iterative way
% Bg : Background, assumed to be known exactly
% mystart = [Intensity, PosX, PosY, SigX, SigY]  : Start vector of intensity, positions and widths
% method : 'COM' uses center of mass and StdDev, 'Int' linear method with intensity weighting, 'IntBg' linear method with intensity and background weighting

function [spots, res, fitted, residual] = FitOneGaussFast(Bg, Meas, method)
if nargin < 3
    method = 'IntBg';
end
EPS=0.001;
Meas2=Meas-Bg;
Mask=Meas2<=EPS;
Meas2(Mask)=EPS;  % to avoid problems with the logarithm
LogI = double(log(Meas2(:)));

% Calculate the right hand side vector
% M * coeff = LogI
% coeff = [e a b c d]
% M = [1 ... 1_n; x1...x_n; x1^2 ... xn^2;y1 ..yn; y1^2...yn^2]
MyXX=xx(Meas2,'corner');
MyXX=double(MyXX(:));  % flatten
MyYY=yy(Meas2,'corner');
MyYY=double(MyYY(:));  % flatten

switch method
    case 'COM'
    % Meas2=Meas-Bg;  % Do not CLIP?
    Meas2 = double(Meas2(:));
    sMeas = sum(Meas2);
    PosX = sum(MyXX .* Meas2)/sMeas;
    PosY = sum(MyYY .* Meas2)/sMeas;
    sX2 = sum(Meas2 .* (MyXX-PosX).^2)/sMeas;
    sY2 = sum(Meas2 .* (MyYY-PosY).^2)/sMeas;
    MyInt = sMeas;
        
    case {'Int', 'IntBg'}
        if strcmp('Int',method)
            weights = Meas2.^2 / (Meas2); 
        else
            weights = Meas2.^2 / (Meas2 + Bg); 
        end
    weights = double(weights(:));
    M = [ones(1,size(LogI,2)); MyXX; MyXX.^2; MyYY; MyYY.^2] .* repmat(weights(:),[1 5])';
    % Calculate the equation system to invert

    sol=pinv(M)' * (LogI .* weights)';

    sX2 = -1/(2*sol(3));
    sY2 = -1/(2*sol(5));
    PosX = sol(2)*sX2;
    PosY = sol(4)*sY2;
    MyInt = exp(sol(1) + PosX^2/(2*sX2) + + PosY^2/(2*sY2));
    otherwise
        Error('unknow fitting method');
end
spots=[MyInt PosX PosY sqrt(sX2) sqrt(sY2)];
