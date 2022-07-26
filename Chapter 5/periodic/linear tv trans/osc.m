clc;  clear all;


x = [.5;.1];
X0 = [.45, .5; ...
      .1, .2];



Tf = 10;
[t_x, xs] = ode45(@(t,y) F(y), [0, Tf], x);
t_x = t_x';
xs = xs';

[t_a, as] = ode45(@(t,a) E(a), [0, Tf], X0(:));
t_a =t_a';
as = as';

options = optimoptions('fsolve','Display','none');
aeq = fsolve(@E, [0; -0.5; .6; 1], options)

figure(1); clf;hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
%axis([-.1, .6, -1.1, 1.5])
xlabel('$V_{\rm d}$','Interpreter','latex')
ylabel('$i_{\rm L}$','Interpreter','latex')

points = makeRectangle(aeq(:, 1));
patch(points(1, :), points(2, :), 'g', 'LineWidth',1, 'FaceAlpha', .3);

points = makeRectangle(as(:, 1));
patch(points(1, :), points(2, :), 'r', 'LineWidth',1, 'FaceAlpha', .3);



[q1, q2] = meshgrid(X0(1, 1):(X0(1, 2) - X0(1, 1))/3:X0(1, 2), ...
                    X0(2, 1):(X0(2, 2) - X0(2, 1))/3:X0(2, 2));
holder = [q1(:), q2(:)]';
for i = 1:size(holder, 2)
    x = holder(:, i);
    [t_x, xs] = ode45(@(t,y) F(y), [0, Tf], x);
    t_x = t_x';
    xs = xs';
    plot(xs(1, :), xs(2, :), 'b', 'LineWidth',2);
    drawnow
end
%options = optimoptions('fsolve','Display','none');
%aeq = fsolve(@E2, [.2; 0.1; .3; .2], options)


function out = F(x)
    C = 1;  %10^(-12);
    L = 1;  %10^(-6);
    R = 1/5;% 200;

    out(1, 1) = 1/C*(- ID(x(1)) + x(2));
    out(2, 1) = 1/L*(-x(1) - R*x(2) + .3);
end

function out = ID(vd)
    out = 803.712*vd^5 - 1086.288 *vd^4 + 551.088*vd^3 - 124.548 *vd^2 + 10.656*vd;
end


function out = E(a)
    x = a(1:2, 1);
    xh = a(3:4, 1);
    out = [d(x, xh); d(xh, x)];
end

function out = E2(a)
    x = a(1:2, 1);
    xh = a(3:4, 1);
    out = [d2(x, xh); d2(x, xh)];
end

function out = d(x, xh)
    C = 1;  %10^(-12);
    L = 1;  %10^(-6);
    R = 1/5;% 200;
    out(1, 1) = 1/C*(- ID(x(1)) + x(2));
    out(2, 1) = 1/L*(-xh(1) - R*x(2) + .3);
end

function out = d2(x, xh)
    C = 1;  %10^(-12);
    L = 1;  %10^(-6);
    R = 1/5;% 200;
    out(1, 1) = 1/C*(ID(x(1)) - xh(2));
    out(2, 1) = 1/L*(x(1) + R*x(2) - .3);
end

function out = makeRectangle(Xin)
    X = Xin(:);
    out = [X(1), X(3), X(3), X(1); ...
           X(2), X(2), X(4), X(4)];
end





