% myplane=hyperplaneString(mydim,mypos,numdims) : Generates the argument string for selecting a hyperplane at position mypos (0-based notation),orthogonal to the dimension mydim 
% This can be used for acces and assignment
% example:
% eval(['hyperplane=myimage' hyperplaneString(2,18,3) ';'])

function s=hyperplaneString(mydim,mypos,numdims)
    s='(';
    for d=1:numdims
        if d== mydim
            s=[s num2str(mypos)];
        else
            s=[s ':'];
        end
        if d<numdims
            s=[s ','];
        end
    end
    s=[s ')'];
end