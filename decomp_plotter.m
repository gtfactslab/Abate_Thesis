


global a b

a = 1;




dx = .2;
X_range = [-2, 2];

[X, Xh] = meshgrid(X_range(1):dx:X_range(2), ...
                   X_range(1):dx:X_range(2));
num = size(X, 1);
               
points = [X(:)'; Xh(:)'];


real = [];
% get diagnol
for x = X_range(1):dx:X_range(2)
    out = F(x);
    real = [real, [x; x; out]];
end

holder = [];
for i = 1:size(points, 2)
    xnow = points(:, i);
    out = d(xnow(1), xnow(2));
    out2 = d2(xnow(1), xnow(2));
    
    holder = [holder, [xnow; out; out2]];    
end
X = reshape(holder(1, :), [num, num]);
Y = reshape(holder(2, :), [num, num]);
Z1 = reshape(holder(3, :), [num, num]);
Z2 = reshape(holder(4, :), [num, num]);

figure(1); clf;
hold on; grid on;
axis([X_range(1), X_range(2), ...
      X_range(1), X_range(2), ...
      3*[-1, 1]])
xlabel('x')
ylabel('xh')
zlabel('d')
view([-.5, 1, 1])
surf(X, Y, Z1)
surf(X, Y, Z2, 'FaceColor', 'r')
plot3(real(1, :),real(2, :),real(3, :), 'r', 'LineWidth', 2);

function out = F(x)
    global a b
    out = a*x^2;
end



function out = d(x, xhat)
    if x <= xhat
        fun_y1 = @(x) F(x);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, x, xhat, options);
        
        out(1, 1) = fun_y1(out_y1);
        
    elseif xhat <= x
        fun_y1 = @(x) - F(x);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, xhat, x, options);
       
        out(1, 1) = - fun_y1(out_y1);
    end 
end

function out = d2(x, xhat)
    global a b
    
    if a*x >= max(0, -a*xhat)
        out = a*x^2;
    elseif a*xhat <= min(0, -a*x)
        out = a*xhat^2;
    elseif a*x <= 0 && 0 <= a*xhat 
        out = 0;
    end 
end