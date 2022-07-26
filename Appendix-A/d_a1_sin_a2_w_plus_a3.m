function out = d_a1_sin_a2_w_plus_a3(w, what, a)
    a1 = a(1);
    a2 = a(2);
    a3 = a(3);

    if a1*a2 >= 0
        out = a1*d_sin(a2*w + a3, a2*what + a3);
    elseif a1*a2 < 0
        out = a1*d_sin(a2*what + a3, a2*w + a3);
    end
end