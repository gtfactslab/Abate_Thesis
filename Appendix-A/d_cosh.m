function out = d_cosh(w, what)
    if w >= max([0, -what])
        out = cosh(w);
    elseif what <= min([0, -w])
        out = cosh(what);
    else
        out = 1;
    end
end