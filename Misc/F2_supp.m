

clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global W
% Bound on disturbance
W = [-3/4, 3/4; ...
     -3/4, 3/4];
 
dx = 1
[X,Y] = meshgrid(-15:dx:15 , -15:dx:15);
Z1 = zeros(size(-15:dx:15, 2));
Z2 = zeros(size(-15:dx:15, 2));

for i = 1:size(Z1, 1)
    for j = 1:size(Z1, 2)
        Z1(i, j) = d(X(i, j), Y(i, j));
        
        Z2(i, j) = q(X(i, j), Y(i, j));
    end
end

figure(1); clf;
hold on; grid on;
surf(X,Y,Z1)
surf(X,Y,Z2, 'FaceColor', 'g')
view([-.7 1 .5])


%%%%%%%%%%%%%%%%%%%%
% without transformation
%%%%%%%%%%%%%%%%%%%%
% dynamics 

function out = F(x)
    c = -1;
    out = -x - c * x^2;
end 


% forward-time decomposition function
function out = q(a, b)
    c = -1;
    
    if c > 0
        if c*a <= -1/2 && -1/2 <= c*b
            out = min([F(a), F(b)]);
        end
    end    
    if c < 0
        if c*a <= -1/2 && -1/2 <= c*b
            out = max([F(a), F(b)]);
        end
    end
    
    if c*b <= -1/2 && -1/2 <= c*a
        out = 1/(4*c);
    end
    
    if -1/2 >= max([c*a, c*b]) 
        out = F(a);
    end
    
    if -1/2 <= min([c*a, c*b])
         out = F(b);
    end
    

end

% forward-time decomposition function
function out = d(y, yhat)
    
    if y <= yhat
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) F(x);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y, yhat, options);
        
        out(1, 1) = fun_y1(out_y1);
        
    elseif y >= yhat
        fun_y1 = @(x) -F(x);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat, y, options);
        
        out(1, 1) = - fun_y1(out_y1);
    else
        out = [inf; inf];
    end
end


