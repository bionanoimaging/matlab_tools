% myplane=hyperplane(img,myplane,mydim,mypos) : Assigns a hyperplane to the data img at position mypos (0-based notation),orthogonal to the dimension mydim 

function s=hyperplaneAsgString(myplane,mydim,mypos)
    s='img(';
    for d=1:ndims(img)
        if d== mydim
            s=[s num2str(mypos)];
        else
            s=[s ':'];
        end
        if d<ndims(img)
            s=[s ','];
        end
    end
    s=[s ')=myplane;'];
end