% Paper:       "Tight Decomposition Functions for Continuous-Time 
%               Mixed-Monotone Systems with Disturbances"
% Publisher:    Control Systems Letters (L-CSS)
% 
% Description: This script generates Figure 1 in the paper. Two tight
%              decomposition functions are constructed for a dynaical 
%              system.  It is shown that the tight decomposition function 
%              proposed in the paper leads to tighter approximations.
%
% Code Author: Matthew Abate
% Date:        6/5/2020

clc; clear all


X0 = [-1, 1; ...
      0, 1]
disc = 20;
holder1 = linspace(X0(1, 1), X0(1, 2), disc);
holder2 = linspace(X0(2, 1), X0(2, 2), disc);
X = [];
for i =1:1:disc
    for j =1:1:disc
        for k =1:1:disc
            X = [X, [holder1(i); holder2(j)]];
        end
    end
end
clear holder1 holder2

k = boundary(X(1, :)', X(2, :)');
X = X(:, k); 
X = X(:);

dt = .002;
T = 0: dt :1;
for t = T
    for j = 1:size(X, 1)/2
        xnow = X([2*j - 1, 2*j], 1);
        xnext = xnow + dt* F(xnow);
        X([2*j - 1, 2*j], 1) = xnext;
    end
end

X = reshape(X, 2, size(X, 1)/2);



figure(1); clf; hold on; grid on; 
axis([-1.5 4.5 -2.5 2.5])
xlabel('$x_1$','Interpreter','latex')
xticks([-1.5, 0, 1.5, 3, 4.5])
ylabel('$x_2$','Interpreter','latex')
yticks([-2.5, 0, 2.5])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')



b = X0(:);
for t = T
    for j = 1:size(X, 1)/2
        b = b + dt* [d2(b(1:2, 1), b(3:4, 1)); d2(b(3:4, 1), b(1:2, 1))];
    end
end
b = reshape(b, 2, 2);
over_approx = [b(1, 1), b(1, 2),b(1, 2),b(1, 1); ...
               b(2, 1), b(2, 1),b(2, 2),b(2, 2)];
patch(over_approx(1, :), over_approx(2, :), 'm', 'FaceAlpha', .2, 'LineWidth', 1.5)


a = X0(:);
for t = T
    for j = 1:size(X, 1)/2
        a = a + dt* [d(a(1:2, 1), a(3:4, 1)); d(a(3:4, 1), a(1:2, 1))];
    end
end
a = reshape(a, 2, 2);
over_approx = [a(1, 1), a(1, 2),a(1, 2),a(1, 1); ...
               a(2, 1), a(2, 1),a(2, 2),a(2, 2)];
patch(over_approx(1, :), over_approx(2, :), 'w', 'LineWidth', 1.5)
patch(over_approx(1, :), over_approx(2, :), 'b', 'FaceAlpha', .2, 'LineWidth', 1.5)

patch(X(1, :), X(2, :), 'w', 'FaceAlpha', 1, 'LineWidth', 1.5)
b = X0;
over_approx = [b(1, 1), b(1, 2),b(1, 2),b(1, 1); ...
               b(2, 1), b(2, 1),b(2, 2),b(2, 2)];
patch(over_approx(1, :), over_approx(2, :), 'r', 'LineWidth', 1.5)

patch(X(1, :), X(2, :), 'g', 'FaceAlpha', .8, 'LineWidth', 1.5)

Leg = legend();
set(Leg,'visible','off')

%matlab2tikz('corbin.tikz', 'width', '5.5cm', 'height', '3.5cm')

function dxdt = F(x)
    dxdt = [abs(x(1) - x(2)); ...
            -x(1)];
end

function out = d(x, xh)
    if x(2) <= x(1) && x(1) <= xh(2)
        out(1, 1) = 0;
    elseif ( 2*x(1) <= 2*x(2)) && (2*x(1) <= xh(2) + x(2))
        out(1, 1) = x(2)- x(1);
    elseif ( 2*x(1) >= xh(2)) && (2*x(1) >= xh(2) + x(2))
        out(1, 1) = x(1) - xh(2);
    end
    out(2, 1) = - xh(1);
end


function out = d2(x, xh)
    if x(1) <= xh(1) && xh(1) <= x(2) && x(2) <= xh(2)
        out(1, 1) = x(2) - xh(1);
        1
    elseif x(2) <= xh(2) && xh(2) <= x(1) && x(1) <= xh(1)
        out(1, 1) = x(1)- xh(1);
        2
    elseif x(2) <= xh(1) && x(1) <= xh(2)
        out(1, 1) = 0;
        3
    elseif xh(1) <= min([x(1), xh(2)]) && max([x(1), xh(2)]) <= x(2)
        out(1, 1) = x(2) - xh(1);
        4
    elseif xh(2) <= min([x(2), xh(1)]) && max([x(2), xh(1)]) <= x(2)
        out(1, 1) = x(1) - xh(2);
        5
        
    elseif xh(1) <= xh(2) && xh(2) <= x(2) && x(2) <= x(1)
        out(1, 1) = x(2) - xh(1);
        6
    elseif xh(2) <= xh(1) && xh(1) <= x(1) && x(1) <= x(2)
        out(1, 1) = x(1)- xh(1);
        7
    end
    out(2, 1) = - xh(1);
end

