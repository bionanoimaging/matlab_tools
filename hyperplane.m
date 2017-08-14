% myplane=hyperplane(img,mydim,mypos) : Extracts a hyperplane from the data img, orthogonal to the dimension mydim and at position mypos (0-based notation)

function myplane=hyperplane(img,mydim,mypos)
    s='myplane=img(';
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
    s=[s ');'];
    eval(s);
end