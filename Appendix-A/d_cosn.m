function out = d_cosn(w, what, a1)

    if mod(a1, 2) == 1
        out = (d_cos(w, what))^(a1);
    elseif mod(a1, 2) == 0
        out = (d_cos2(w, what))^(a1/2);
    elseif a1 >= 0 && -pi/2 <= min([w, what]) && min([w, what]) < pi/2
        out = (d_cos(w, what))^(a1);
    elseif a1 <= 0 && -pi/2 < min([w, what]) && min([w, what]) < pi/2
        out = (d_cos(what, w))^(a1);
    else
        error('Error: d_sinn requires integer n, or w, what in [0, pi]');
    end
end