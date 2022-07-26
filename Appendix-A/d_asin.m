function out = d_asin(w, what)
    if -1 <= min([w, what]) && max([w, what]) <= 1
        out = asin(w);
    else
        error('Error: d_asin domain is [-1, 1]');
    end
end