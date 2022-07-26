clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global A B

A = [0, -1; ...
     1, 0];


X0 = [2, 2.5; ...
      1, 1.5];

B = [0; 0];
W = [0, 1];

T_sim = 1;
dx = .1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XE, XE_Traj] = PhiE(T_sim, [X0(:, 1); X0(:, 2)], W(:), A, B);


X_initial = rect4plot(X0);

Reachable_Set = Reach(T_sim, X0, W);
k = boundary(Reachable_Set(1, :)', Reachable_Set(2, :)');
Reachable_Set = Reachable_Set(:, k);

% plot MM estimate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: Just Reachable Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
axis([-1, 3, 0, 4])
xticks(-1:1:3);
yticks(0:2:4);

% plot initial set
X0_plot = rect4plot(X0);
patch(X0_plot(1, :), X0_plot(2, :), 'r', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');

XE_plot = rect4plot([XE(1), XE(3); XE(2), XE(4)]);
patch(XE_plot(1, :), XE_plot(2, :), 'b', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');


y0 = X0(:);
dt = .01;
T_vec = 0:dt:T_sim ;
y_now = y0;
for i = 1:size(T_vec, 2)
    t_now = T_vec(i);
    T_now_inv = expm(-A*t_now);
    BM = T_now_inv*B;
    BMp = BM.*(BM >= 0);
    BMm = BM - BMp;
    y_now = y_now + dt*[BMp, BMm; BMm, BMp]*W(:);
end
yt = y_now;
Tt = expm(A*t_now);
XE_plot = Tt*rect4plot(reshape(yt, [2, 2]));

patch(XE_plot(1, :), XE_plot(2, :), 'r', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');

ApM =  Tt.*(Tt >= 0);
AmM = Tt - ApM;

AM = [ApM, AmM ; AmM, ApM];
Y0 = AM*yt;


XE_plot = rect4plot([Y0(1), Y0(3); Y0(2), Y0(4)]);
patch(XE_plot(1, :), XE_plot(2, :), 'r', 'FaceColor', 'none', ...
                            'LineStyle', '--', ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');



% plot true reachable set
holder = Reachable_Set(:, 1);
for i = 2:size(Reachable_Set, 2)
    if norm(holder(:, end) - Reachable_Set(:, i)) >= .01
        holder(:, end+1) =  Reachable_Set(:, i);
    else
    end
end

patch(holder(1, :), holder(2, :), 'w')
patch(holder(1, :), holder(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off')


% plot points
points = [X0, XE(1:2, 1), XE(3:4, 1), Y0(1:2, 1), Y0(3:4, 1)]
scatter([points(1, :)], points(2, :), 'black', 'filled')

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

matlab2tikz('linear_reach2.tikz', 'width', '6cm', 'height', '4cm')

function out = rect4plot(X)
    out = [X(1, 1), X(1, 2), X(1, 2), X(1, 1); ...
           X(2, 1), X(2, 1), X(2, 2), X(2, 2)];
end

function out = Decomposition(X, U, A, B)
    Ap =  (A - diag(diag(A))).*(A >= 0) + diag(diag(A));
    Am = A - Ap;
    Bp =  B.*(B >= 0);
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
    Xflip = [X(end/2+1:end, 1); X(1:end/2, 1)];
    Uflip = [U(end/2+1:end, 1); U(1:end/2, 1)];
    
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
    dt = .001;
    T_sim = 0:dt:T;
    U_sim = U(1):(U(2) - U(1))/2:U(2);
    
    [X1, X2] = meshgrid(X0(1, 1):(X0(1, 2) - X0(1, 1))/5: X0(1, 2), ... 
                        X0(2, 1):(X0(2, 2) - X0(2, 1))/5: X0(2, 2))

    initial_set = [X1(:), X2(:)]';
    k = boundary(initial_set(1, :)', initial_set(2, :)');
    initial_set = initial_set(:, k);

    Reachable_Set = initial_set;
    for i = 1:size(T_sim, 2) -1
        i
        X_Next = [];
        for j = 1:size(Reachable_Set, 2)
            x_now = Reachable_Set(:, j);
            for u = U_sim
                x_next = x_now + dt*Dynamics(x_now, u);
                X_Next = [X_Next, x_next];
            end
        end
        k = convhull(X_Next(1, :)', X_Next(2, :)' );
        Reachable_Set = X_Next(:, k);
    end
end

function dxdt = Dynamics(x, u)
    global A B
    dxdt = A*x + B*u;
end






