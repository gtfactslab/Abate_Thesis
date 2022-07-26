
% Title: Decomposition Functions for Interconnected Mixed Monotone Systems
% Submitted to IEEE Control Systems Letters (L-CSS), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 11/14/2021
% Description:  This script generates Figure 3.
%               A tight decomposition function for the kinematic unicycle 
%               model is formed by applying Theorem 2, and a 
%               hyperrectangular approximation of the systems reachable set 
%               is computed by applying Proposition 1.


clc; clear all;

X0 = [-1, 1; ...
      -1, 1; ...
      pi/3, pi/2];

U = [-.1, .1];
du1 = .05;
du2 = .05;

dx = .05;
d0 = .01;

l = 1;

% set up
X0_set = zeros(3, 4*length(X0(1, 1):dx:X0(1, 2)) - 4);
for i = X0(1, :)
    for j = X0(2, 1):dx:X0(2, 2) 
        holder(:, l) = [i; j];
        l=l+1;
    end
end
for i = X0(1, 1):dx:X0(1, 2) 
    for j = X0(2, :) 
        holder(:, l) = [i; j];
        l=l+1;
    end
end
k = convhull(holder(1, :)', holder(2, :)');
holder = holder(:, k(1:end-1));
l=1;
for i = 1:size(holder, 2)
    for k = X0(3, 1):d0:X0(3, 2) 
        X0_set(:, l) = [holder(:, i);k];
        l = l+1;
    end
end

l = 1;
for i = U(1, 1):du1:U(1, 2)
        U_set(:, l) = [i];
        l = l+1;
end


thing = [];
l = 1;
T_sim = 6;
for i = 1:length(X0_set)
    for j = 1:size(U_set, 2)
        thing(:, l) = sim(T_sim, X0_set(:, i), U_set(:, j));
        l = l + 1;
    end
end

%
figure(1); clf;
hold on; grid on;
Leg = legend();
set(Leg,'visible','off');
axis([-5, 7, -6, 15])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$p_{x}$','Interpreter','latex')
ylabel('$p_{y}$','Interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
over_approx = sim_E(T_sim, X0(:), U, @d_tight);
over_approx = reshape(over_approx, [3, 2]);

l = 1;
for i = over_approx(1, :)
    for j = over_approx(2, :)
        corn(:, l) = [i, j];
        l = l+1;
    end
end
k = convhull(corn(1, :)', corn(2, :)');
patch(corn(1, k), corn(2, k), 'w', 'FaceAlpha', 1);
patch(corn(1, k), corn(2, k), 'g', 'FaceAlpha', .1);

k = boundary(thing(1, :)', thing(2, :)', .3);
patch(thing(1, k(1:1:end)), thing(2, k(1:1:end)), 'g');


patch([-1, 1, 1, -1], [-1, -1, 1, 1], 'w', 'FaceAlpha', .8);
patch([-1, 1, 1, -1], [-1, -1, 1, 1], 'r', 'FaceAlpha', .9);
drawnow;

matlab2tikz('unicycle_3.tikz', 'width', '6cm', 'height', '4cm')

function xdot = F(x, w)
    px_dot = cos(x(3));
    py_dot = sin(x(3));
    theta_dot = w;

    xdot = [px_dot;py_dot;theta_dot];
end


function out = sim(T, x0, u)
    X = x0;
    dt = T/1000;
    timevec = 0:dt:T;
    for t = 1:length(timevec)
        X = X + dt*F(X, u);
    end
    out = X;
end

function out = sim_E(T, x0, U, d)
    X = x0;
    dt = T/1000;
    timevec = 0:dt:T;
    for t = 1:length(timevec)
        X = X + dt*E(X(1:3), U(:, 1), X(4:6), U(:, 2), d);
    end
    out = X;
end

function out = E(x, u, xh, uh, d)
    out = [ d(x, u, xh, uh); ...
            d(xh, uh, x, u) ] ;
end

function out = d_norm(x, w, xh, wh)
    out = [ cos(x(3)) + x(3) - xh(3); ...
            sin(x(3)) + x(3) - xh(3); ...
            w];
end

function out = d(x, w, xh, wh)
    out = [ d1(x, w, xh, wh); ...
            d2(x, w, xh, wh); ...
            w];
end



function out = d1(x, w, xh, wh)
    out = 10000;
    if all([x; w] <= [xh; wh])
        out = max(cos(x(3)) + x(3) - xh(3), cos(xh(3)) + x(3) - xh(3));
    else
        out = min(cos(x(3)) + x(3) - xh(3), cos(xh(3)) + x(3) - xh(3));
    end
end

function out = d2(x, w, xh, wh)
    out = 10000;
    if all([x; w] <= [xh; wh])
        out = max(cos(x(3)) + x(3) - xh(3), ...
                  cos(xh(3)) + x(3) - xh(3));
    else
        out = min(sin(x(3)) + x(3) - xh(3), ...
                  sin(xh(3)) + x(3) - xh(3));
    end
end

function out = d_best(x, w, xh, wh)
    out = [ max(min(1, d1(x, w, xh, wh)), -1); ...
            max(min(1, d2(x, w, xh, wh)), -1); ...
            w];
end

function out = d_tight(x, w, xh, wh)
    out = [ d_cos(x(3), xh(3)); ...
            d_sin(x(3), xh(3)); ...
            w];
end




