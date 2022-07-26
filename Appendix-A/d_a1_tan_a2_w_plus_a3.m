function out = d_a1_tan_a2_w_plus_a3(w, what, a)
    a1 = a(1);
    a2 = a(2);
    a3 = a(3);
    Domain = [min([(-pi - 2*a3)/(2*a2), (pi - 2*a3)/(2*a2)]), ...
              max([(-pi - 2*a3)/(2*a2), (pi - 2*a3)/(2*a2)])];

    if Domain(1) < min([w, what]) && max([w, what]) < Domain(2)
        if a1*a2 >= 0
            out = a1*tan(a2*w + a3);
        elseif a1*a2 < 0
            out = a1*tan(a2*what + a3);
        end
    else
        error('Error: d_a1_tan_a2_w_plus_a3 domain is ((-pi - 2 a3)/(2 a2), (pi - 2 a3)/(2 a2))');
    end
end