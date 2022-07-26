function out = d_asin(w, what)
    if -1 <= min([w, what]) && max([w, what]) <= 1
        out = acos(what);
    else
        error('Error: d_acos domain is [-1, 1]');
    end
end