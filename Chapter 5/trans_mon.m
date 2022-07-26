
% Title: Improving the Fidelity of Mixed-Monotone Reachable Set 
%        Approximations via State Transformations
% Conference: American Controls Conference (ACC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 3/18/2021
% Description:  This script generates Figure 3.
%               A linear transformation is found to transform a given
%               system to a monotone system.  This allows for the
%               computation of the tightest parallelogram containing the
%               reachable set of the system via Theorem 1.  A second
%               approximation is also computed via a different
%               transformation, and this allows for reduced conservatism in
%               the approximation.


clc; clear all;

% simulation parameters
dt = .01;
Sim_Time = 1;

global T1 T2 W
T1 = [ 1, 1/sqrt(2) ;...
      0,  1/sqrt(2) ];
  
T2 = [  1, 2/sqrt(3) ;...
       0, 1/sqrt(3) ];
  
W = [-1, 2]

syms x1 x2 w y1 y2 a b

T = [1, a; 0, 1];
simplify(dydt([x1; x2], w, T))
%Jy = [ diff(dydt([y1; y2], w, T), y1), diff(dydt([y1; y2], w, T), y2)];
%simplify(Jy)

%


x0 = [1; 0];
x_E1 = [x0; x0];
x_E2 = [inv(T1)*x0; inv(T1)*x0];
x_E3 = [inv(T2)*x0; inv(T2)*x0];


holder = x0;
for t = 0:dt:Sim_Time
    holder2 = [];
    for j = 1:1:size(holder, 2)
        xnow = holder(:, j);
        for w = W(1):(W(2)-W(1))/4:W(2)
            xnext = xnow + dt*F(xnow, w);
            holder2 = [holder2, xnext];
        end
    end
    if t > 2*dt
        k = boundary(holder2(1, :)', holder2(2, :)');
        holder = holder2(:, k);
    else
        holder = holder2;
    end

    x_E1 = x_E1 + dt*E(x_E1(1:2, 1), x_E1(3:4, 1), @d1); % rect    
    x_E2 = x_E2 + dt*E(x_E2(1:2, 1), x_E2(3:4, 1), @d2); % parallelogram 1
    x_E3 = x_E3 + dt*E(x_E3(1:2, 1), x_E3(3:4, 1), @d3); % parallelogram 1
end
ReachT = holder;
clear k holder xnow xnext

%
figure(1); clf;
hold on; grid on;
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1, 6, -.25, 1.75]);
xticks(0:2:6)
yticks(0:.5:1.5)


scatter(x0(1), x0(2), 80, 'r', 'filled')

ETT_Approx = makeRectangle([x_E1(1:2), x_E1(3:4)]);
patch(ETT_Approx(1, :), ETT_Approx(2, :), 'b', 'FaceAlpha', .1, 'LineWidth', 1.25)

ET_Approx = T1*makeRectangle([x_E2(1:2), x_E2(3:4)]);
patch(ET_Approx(1, :), ET_Approx(2, :), 'w', 'FaceAlpha', .2, 'LineWidth', 1.25)
patch(ET_Approx(1, :), ET_Approx(2, :), 'r', 'FaceAlpha', .1, 'LineWidth', 1.25)


ET_Approx = T2*makeRectangle([x_E3(1:2), x_E3(3:4)]);
patch(ET_Approx(1, :), ET_Approx(2, :), 'w', 'FaceAlpha', .2, 'LineWidth', 1.25)
patch(ET_Approx(1, :), ET_Approx(2, :), 'r', 'FaceAlpha', .1, 'LineWidth', 1.25)


%plot tru reachable set, from exhaustive simulation

holder = ReachT(:, 1);
for i = 2:size(ReachT, 2)
    if norm(holder(:, end) - ReachT(:, i)) >= .11
        holder(:, end + 1) = ReachT(:, i);
    end
end
patch(holder(1, :), holder(2, :), 'w', 'FaceAlpha', 1, 'LineWidth', 1.25)
patch(holder(1, :), holder(2, :), 'g', 'FaceAlpha', 1, 'LineWidth', 1.25)

Leg = legend();
set(Leg,'visible','off')
grid on;
ax.Layer = 'top';

drawnow

% 
% % true reachable set
% shp1 = alphaShape(ReachT(1, 1:disc:end)',ReachT(2, 1:disc:end)');
% area_reach(1) = area(shp1);
% 
% area_reach(2) = prod(max(ET_Approx')' - min(ET_Approx')');
% 
% % approximation from T1
% shp2 = alphaShape(ET_Approx(1, :)', ET_Approx(2, :)');
% shp2.Alpha = 2.5;
% area_reach(3) = area(shp2);
% 
% % true reachable set
% shp3 = alphaShape(ETT_Approx(1, :)', ETT_Approx(2, :)');
% area_reach(4) = area(shp3);
% 
% area_reach
% 
% 




matlab2tikz('trans_mon.tikz', 'width', '6cm', 'height', '4cm')

function out = F(x, w)
    out = [ x(1) - 2*x(2) + x(2)^3 + w; ...
            (x(1) - 2*x(2))];
end

function out = dydt(x, w, T1)
    %global T1
    out = inv(T1)*F(T1*x, w);
end

function out = dydt2(x, w, T2)
    %global T1
    out = inv(T2)*F(T2*x, w);
end

function out = d1(x, w, xh, wh)
    if x(1) <= xh(1) && x(2) <= xh(2) && w <= wh
            holder = [- 2*x(2) + x(2)^3, - 2*xh(2) + xh(2)^3];
            if x(2) <= -sqrt(2/3) && -sqrt(2/3) <= xh(2)
                holder = [holder, - 2*(-sqrt(2/3)) + (-sqrt(2/3))^3];
            end
            if x(2) <= sqrt(2/3) && sqrt(2/3) <= xh(2)
                holder = [holder, - 2*(sqrt(2/3)) + (sqrt(2/3))^3];
            end
            out(1, 1) = x(1) + w + min(holder);
            out(2, 1) = x(1) - 2*x(2);

    elseif xh(1) <= x(1) && xh(2) <= x(2) && wh <= w
        holder = [- 2*x(2) + x(2)^3, - 2*xh(2) + xh(2)^3];
            if xh(2) <= -sqrt(2/3) && -sqrt(2/3) <= x(2)
                holder = [holder, - 2*(-sqrt(2/3)) + (-sqrt(2/3))^3];
            end
            if xh(2) <= sqrt(2/3) && sqrt(2/3) <= x(2)
                holder = [holder, - 2*(sqrt(2/3)) + (sqrt(2/3))^3];
            end
            out(1, 1) = x(1) + w + max(holder);
            out(2, 1) = x(1) - 2*x(2);
    else
        out = [inf; inf];
    end
end


function out = d2(x, w, xhat, what)
    global T1
%     out = [ x(2)^3 + w; ...
%             x(1) - x(2)];
      out = inv(T1)*F(T1*x, w);
end

function out = d3(x, w, xhat, what)
%     out = [ x(2)^3 - x(1) + w; ...
%             x(1) ];
     global T2
      out = inv(T2)*F(T2*x, w);
end

function out = E(x, xhat, d)
    global W
    out(1:2, 1) = d(x, W(1), xhat, W(2));
    out(3:4, 1) = d(xhat, W(2), x, W(1));
end

function out = makeRectangle(X0)
    out = [X0(1, 1), X0(1, 2), X0(1, 2), X0(1, 1); ...
           X0(2, 1), X0(2, 1), X0(2, 2), X0(2, 2)];
end






