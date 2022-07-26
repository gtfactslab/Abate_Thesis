clc;  clear all;


x = [.5;.1];
X0 = [.45, .5; ...
      .1, .1];
W = [.3, .3];



Tf = 3;
[t_x, xs] = ode45(@(t,y) F(y, W(1)), [0, Tf], x);
t_x = t_x';
xs = xs';

[t_a, as] = ode45(@(t,a) E(a, W), [0, Tf], X0(:));
t_a =t_a';
as = as';

options = optimoptions('fsolve','Display','none');
aeq = fsolve(@(a) E(a, W), [0; 0; 2; 2], options);


figure(1); clf;hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
%axis([-.1, .6, -1.1, 1.5])
xlabel('$V_{\rm d}$','Interpreter','latex')
ylabel('$i_{\rm L}$','Interpreter','latex')

C = 1;  %10^(-12);
L = 1;  %10^(-6);
R = 1/5;% 200;
[t_a, as] = ode45(@(t,a) Et(t, a), [0, Tf], X0(:));
t_a =t_a';
as = as';
for i = fliplr(1:size(t_a, 2))
    t = t_a(i);
    points = [1, 0; 0, exp(-R/L*t)]*makeRectangle(as(:, i));
patch(points(1, :), points(2, :), 'r', 'LineWidth',1, 'FaceAlpha', .3);
end


points = makeRectangle(aeq(:, 1));
patch(points(1, :), points(2, :), 'g', 'LineWidth',1, 'FaceAlpha', .3);



points = makeRectangle(as(:, 1));
patch(points(1, :), points(2, :), 'r', 'LineWidth',1, 'FaceAlpha', .3);



plot(xs(1, :), xs(2, :), 'b', 'LineWidth',2);
drawnow
%%
x = sym('x', [2, 1]);
syms w
F2(x, w)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = F(x, w)
    C = 1;  %10^(-12);
    L = 1;  %10^(-6);
    R = 1/5;% 200;
    out(1, 1) = 1/C*(- ID(x(1)) + x(2));
    out(2, 1) = 1/L*(-x(1) - R*x(2) + w);
end

function out = F2(x, w)
    T = [1, 0; ...
         .2, 1];
    out = inv(T)*F(T*x, w);
end


function out = E(a, W)
    x = a(1:2, 1);
    xh = a(3:4, 1);
    out = [d(x, W(:, 1), xh, W(:, 2)); ...
           d(xh,W(:, 2),  x, W(:, 1))];
end

function out = d(x, w, xh, wh)
    C = 1;  %10^(-12);
    L = 1;  %10^(-6);
    R = 1/5;% 200;
    out(1, 1) = 1/C*(- ID(x(1)) + x(2));
    out(2, 1) = 1/L*(-xh(1) - R*x(2) + w);
end



function out = makeRectangle(Xin)
    X = Xin(:);
    out = [X(1), X(3), X(3), X(1), X(1); ...
           X(2), X(2), X(4), X(4), X(2)];
end

function out = Et(t, a)
    x = a(1:2, 1);
    xh = a(3:4, 1);
    out = [dt(t, x, xh); dt(t, xh, x)];
end

function out = dt(t, x, xh)
    C = 1;  %10^(-12);
    L = 1;  %10^(-6);
    R = 1/5;% 200;
    out(1, 1) = 1/C*(- ID(x(1)) + exp((R/L)*t)*x(2));
    out(2, 1) = max(-(1/L)*exp((R/L)*t), 0) * x(1) + ...
                min(-(1/L)*exp((R/L)*t), 0) * xh(1) + ...
                .3*(1/L)*exp((R/L)*t);
end


function out = ID(vd)
    out = 803.712*vd^5 - 1086.288 *vd^4 + 551.088*vd^3 - 124.548 *vd^2 + 10.656*vd;
end