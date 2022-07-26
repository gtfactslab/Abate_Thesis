clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global A B
A = [ -.5,  -1; ...
       1, -1];
B = [2; 0];
W = [1, 1];

T_sim = 2;
dx = .1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
ax = gca;

axis([-2, 4.5, -2, 4.5])
%xticks([-1, 0, 1, 2, 3, 4])
%yticks([-1, 0, 1, 2, 3, 4])
xlabel('$x$','Interpreter','latex')
ylabel('$\widehat{x}$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
% Plot Diagnol of Embedding Space
plot([-2, 6], [-2, 6], 'k', 'LineWidth', 1.25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0 = [-1, 1; -1, 1];
[XE, XE_Traj] = PhiE(T_sim, X0(:), W(:), A, B);
XE_plot = rect4plot([XE(1), XE(3); XE(2), XE(4)]);
xu_T = XE_plot(:, 1);
xo_T = XE_plot(:, 2);

xu = XE_Traj(1:2, :);
xo = XE_Traj(3:4, :);
% Plot trajectory of embedding system that yeilds approximation    
plot(xu(1, :), xo(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 1), X0(1, 2),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xu_T(1, 1), xo_T(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');

% Plot trajectory of embedding system that yeilds approximation    
plot(xo(1, :), xu(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 2), X0(1, 1),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xo_T(1, 1), xu_T(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0 = [-.5, .5; -.5, .5];
[XE, XE_Traj] = PhiE(T_sim, X0(:), W(:), A, B);
XE_plot = rect4plot([XE(1), XE(3); XE(2), XE(4)]);
xu_T = XE_plot(:, 1);
xo_T = XE_plot(:, 2);

xu = XE_Traj(1:2, :);
xo = XE_Traj(3:4, :);
% Plot trajectory of embedding system that yeilds approximation    
plot(xu(1, :), xo(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 1), X0(1, 2),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xu_T(1, 1), xo_T(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');

% Plot trajectory of embedding system that yeilds approximation    
plot(xo(1, :), xu(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 2), X0(1, 1),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xo_T(1, 1), xu_T(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0 = [0, 0; 0, 0];
[XE, XE_Traj] = PhiE(T_sim, X0(:), W(:), A, B);
XE_plot = rect4plot([XE(1), XE(3); XE(2), XE(4)]);
xu_T = XE_plot(:, 1);
xo_T = XE_plot(:, 2);

xu = XE_Traj(1:2, :);
xo = XE_Traj(3:4, :);
% Plot trajectory of embedding system that yeilds approximation    
plot(xu(1, :), xo(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 1), X0(1, 2),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xu_T(1, 1), xo_T(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');


Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

matlab2tikz('symmetry.tikz', 'width', '6cm', 'height', '4cm')

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
    dt = .01;
    T_sim = 0:dt:T;
    U_sim = U(1):.2:U(2);
    
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
        k = convhull(X_Next(1, :)', X_Next(2, :)' );
        Reachable_Set = X_Next(:, k);
    end
end

function dxdt = Dynamics(x, u)
    global A B
    dxdt = A*x + B*u;
end






