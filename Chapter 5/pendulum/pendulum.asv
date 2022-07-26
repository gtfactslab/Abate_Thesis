
clc; clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



kp = 3
kd = 4
global A B
A = [0, 1; ...
     -kp, -kd];

B = [0 ; 1];

X0 = .5*[-1, 1; ...
         -1, 1];

W = [-.5, .5];

T_sim = 7;
dx = .1;




%
%%%%%%%%%%%%%%%%%%%%%

[V, D] = jordan(A)

BD = inv(V)*B;


Dp = (D - diag(diag(D))).*((D - diag(diag(D)))>= 0) + diag(diag(D));
Dm = D - Dp;

DT = [Dp, Dm; Dm Dp];


BDp = (BD >= 0).*BD;
BDm = BD - BDp;
BT = [BDp, BDm; BDm BDp]

aeq = -inv(DT)*BT*W(:);

points = rect4plot([aeq(1:2), aeq(3:4)]);

Xinv = V*points;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XE, XE_Traj] = PhiE(T_sim, X0(:), W(:), A, B);

[a, b] = meshgrid(X0(1, 1):dx:X0(1, 2), ...
                  X0(2, 1):dx:X0(2, 2));

X_initial = [a(:)'; b(:)']; 
clear a b
Reachable_Set = Reach(T_sim, X_initial, W);
k = boundary(Reachable_Set(1, :)', Reachable_Set(2, :)');
Reachable_Set = Reachable_Set(:, k);

X0_plot = rect4plot(X0);

% plot MM estimate
[XT1, X_Traj1] = Phi(T_sim, X0(:, 1), W(:, 1));
[XT2, X_Traj2] = Phi(T_sim, X0(:, 2), W(:, 2));
XT_plot = rect4plot([XT1, XT2]);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: Just Reachable Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
axis([-.5, .5, -.5, .5])
xticks(-.5:.5:.5);
yticks(-.5:.5:.5);

% plot initial set

patch(Xinv(1, :), Xinv(2, :), 'b', ...
                            'FaceAlpha', .2, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off')

Reachable_Set = parse(Reachable_Set, .005);
patch(Reachable_Set(1, :), Reachable_Set(2, :), 'w')
patch(Reachable_Set(1, :), Reachable_Set(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off')



Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';
matlab2tikz('pendulum.tikz', 'width', '6cm', 'height', '4cm')


function out = rect4plot(X)
    out = [X(1, 1), X(1, 2), X(1, 2), X(1, 1); ...
           X(2, 1), X(2, 1), X(2, 2), X(2, 2)];
end

function out = Decomposition(X, U, A, B)
    Ap = [A(1, 1), A(1, 2).*(A(1, 2) >=0); ...
          A(2, 1 ).*(A(2, 1) >=0), A(2, 2)];
    Bp = B.*(B >= 0);
    Am = A - Ap;
    Bm = B - Bp;
   
    out = [Ap, Am]*X + [Bp, Bm]*U;
end

function [XT, XE_Traj] = PhiE(T, X, U, A, B)
    dt = .01;
    T_sim = 0:dt:T ;
    Xt = X;
    for i = 1:size(T_sim, 2)
        x_now = Xt(:, end);
        x_next = x_now + dt*Embedding(x_now, U, A, B);
        Xt = [Xt, x_next];
    end
    XE_Traj = Xt;
    XT = Xt(:, end);
end

function out = Embedding(X, U, A, B)
    Xflip = [X(3:4, 1); X(1:2, 1)];
    Uflip = [U(2, 1); U(1, 1)];
    
    out = [Decomposition(X, U, A, B); 
           Decomposition(Xflip, Uflip, A, B)];
end

function [XT, XE_Traj] = Phi(T, x, u)
    dt = .01;
    T_sim = 0:dt:T ;
    Xt = x;
    for i = 1:size(T_sim, 2)
        x_now = Xt(:, end);
        x_next = x_now + dt*Dynamics(x_now, u);
        Xt = [Xt, x_next];
    end
    XE_Traj = Xt;
    XT = Xt(:, end);
end

function Reachable_Set = Reach(T, X0, U)
    dt = .005;
    T_sim = 0:dt:T;
    U_sim = U(1):.5:U(2);
    
    Reachable_Set = X0;
    for i = 1:size(T_sim, 2)
        i
        X_Next = [];
        for j = 1:size(Reachable_Set, 2)
            x_now = Reachable_Set(:, j);
            for u = U_sim
                x_next = x_now + dt*Dynamics(x_now, u);
                X_Next = [X_Next, x_next];
            end
        end
        k = boundary(X_Next(1, :)', X_Next(2, :)', 0);
        % parse
        Reachable_Set = parse(X_Next(:, k), .002);
    end
end

function out = parse(in, del)
    holder = in;
    holder5 = in(:, 1);
    for j = 2:size(holder, 2)
        if norm(holder5(:, end) - holder(:, j)) >= del
            holder5(:, end + 1) = holder(:, j);
        end
    end
    out = holder5;
end

function dxdt = Dynamics(x, u)
global A B
    dxdt = A*x+B*u;
end






