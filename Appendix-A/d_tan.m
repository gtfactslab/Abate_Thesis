function out = d_tan(w, what)
    if -pi/2 < min([w, what]) && max([w, what]) < pi/2
        out = tan(w);
    else
        error('Error: d_tan domain is (-pi/2, pi/2)');
    end
end