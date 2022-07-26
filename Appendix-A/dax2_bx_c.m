

% x_dot = a1 x^2 + a2 x + a3
function out = dax2_bx_c(x, xh, q)
    a = q(1);
    b = q(2);
    c = q(3);

    if a >= 0 && x <= -b/(2*a) && -b/(2*a) <= xh
        out = (-b^2 +4*a*c)/(4*a);
    elseif a >= 0 && -b/(2*a) <= x && x <= xh
        out = a*x^2 + b*x + c;
    elseif a >= 0 && x <= xh && xh <= -b/(2*a)
        out = a*xh^2 + b*xh + c;
    elseif a >= 0 && xh <= x
        out = max([a*x^2 + b*x + c, a*xh^2 + b*xh + c]);
    
    elseif a <= 0 && xh <= -b/(2*a) && -b/(2*a) <= x
        out = (-b^2 +4*a*c)/(4*a);
    elseif a <= 0 && -b/(2*a) <= xh && xh <= x
        out = a*xh^2 + b*xh + c;
    elseif a <= 0 && xh <= x && x <= -b/(2*a)
        out = a*x^2 + b*x + c;
    elseif a <= 0 && x <= xh
        out = min([a*x^2 + b*x + c, a*xh^2 + b*xh + c]);
    else
        out = nan;
    end
end