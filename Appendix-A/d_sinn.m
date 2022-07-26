function out = d_sinn(w, what, a1)

    if mod(a1, 2) == 1
        out = (d_sin(w, what))^(a1);
    elseif mod(a1, 2) == 0
        out = (d_sin2(w, what))^(a1/2);

    elseif a1 >= 0 && 0 <= min([w, what]) && min([w, what]) <= pi
        out = (d_sin(w, what))^(a1);
    elseif a1 <= 0 && 0 < min([w, what]) && min([w, what]) <= pi
        out = (d_sin(what, w))^(a1);
    else
        error('Error: d_sinn requires integer n, or w, what in [0, pi]');
    end
end