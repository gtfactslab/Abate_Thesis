clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global A B
A = [ -2,  -1; ...
       1, -1];
B = [4; 0];

X0 = .5*[-1, 1; ...
         -1, 1];
W = [1, 2];

T_sim = 1;
dx = .1;


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: Just Reachable Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
axis([-1, 4, -1, 3])
xticks(-1:4);
yticks(-1:3);

% plot initial set
patch(X0_plot(1, :), X0_plot(2, :), 'r');
patch(Reachable_Set(1, :), Reachable_Set(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off')

XE_plot = rect4plot([XE(1), XE(3); XE(2), XE(4)]);
patch(XE_plot(1, :), XE_plot(2, :), 'b', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');

patch(Reachable_Set(1, :), Reachable_Set(2, :), 'w')
patch(Reachable_Set(1, :), Reachable_Set(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off')

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%matlab2tikz('linear_reach1.tikz', 'width', '6cm', 'height', '4cm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: Reachable set but not monotone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); clf; hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
axis([-1, 4, -1, 3])
xticks(-1:4);
yticks(-1:3);

% plot initial set
patch(X0_plot(1, :), X0_plot(2, :), 'r');



scatter(X0(1, 1), X0(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(X0(1, 2), X0(2, 2), 'k', 'filled', 'HandleVisibility', 'off');
scatter(XT1(1, end), XT1(2, end), 'k', 'filled', 'HandleVisibility', 'off');
scatter(XT2(1, end), XT2(2, end), 'k', 'filled', 'HandleVisibility', 'off');

patch(Reachable_Set(1, :), Reachable_Set(2, :), 'w')
patch(Reachable_Set(1, :), Reachable_Set(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off')
patch(XT_plot(1, :), XT_plot(2, :), 'b', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');

plot(X_Traj1(1, :), X_Traj1(2, :), 'b', 'LineWidth', 2);
plot(X_Traj2(1, :), X_Traj2(2, :), 'b', 'LineWidth', 2);
Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%matlab2tikz('linear_reach2.tikz', 'width', '6cm', 'height', '4cm')



%%%%%%%%%%%%%%%%%%%%
figure(3); clf; hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
axis([-1, 4, -1, 3])
xticks(-1:4);
yticks(-1:3);

% plot initial set
patch(X0_plot(1, :), X0_plot(2, :), 'r');

% plot MM estimate
patch(XE_plot(1, :), XE_plot(2, :), 'b', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');
scatter(X0(1, 1), X0(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(X0(1, 2), X0(2, 2), 'k', 'filled', 'HandleVisibility', 'off');
scatter(XE(1, end), XE(2, end), 'k', 'filled', 'HandleVisibility', 'off');
scatter(XE(3, end), XE(4, end), 'k', 'filled', 'HandleVisibility', 'off');

patch(Reachable_Set(1, :), Reachable_Set(2, :), 'w')
patch(Reachable_Set(1, :), Reachable_Set(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off')

plot(XE_Traj(1, :), XE_Traj(2, :), 'b', 'LineWidth', 2);
plot(XE_Traj(3, :), XE_Traj(4, :), 'b', 'LineWidth', 2);

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%matlab2tikz('linear_reach3.tikz', 'width', '6cm', 'height', '4cm')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
hold on; grid on;
ax = gca;
axis([-1, 4, -1, 4])
xticks([-1, 0, 1, 2, 3, 4])
yticks([-1, 0, 1, 2, 3, 4])
xlabel('$x_1$','Interpreter','latex')
ylabel('$\widehat{x}_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Cone Coresponding to X0
patch([X0(1, 1), X0(1, 1), X0(1, 2)], ...
      [X0(1, 2), X0(1, 1), X0(1, 2)], 'r', ...
            'LineWidth', 1, ...
            'FaceAlpha', .7, ...
            'HandleVisibility', 'off');

XE_plot = rect4plot([XE(1), XE(3); XE(2), XE(4)]);
xu_T = XE_plot(:, 1);
xo_T = XE_plot(:, 2);
% Plot Cone Coresponding to Phi^e
patch([xu_T(1, 1), xu_T(1, 1), xo_T(1, 1)], ...
      [xo_T(1, 1), xu_T(1, 1), xo_T(1, 1)], 'g', ...
            'FaceAlpha', .2, ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');

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

% Plot Diagnol of Embedding Space
plot([-1, 4], [-1, 4], 'k', 'LineWidth', 1.25)

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%matlab2tikz('linear_reach4.tikz', 'width', '6cm', 'height', '4cm')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); clf;
hold on; grid on;
ax = gca;
axis([-1, 3, -1, 3])
xticks(-1:3.5)
yticks(-1:3.5)
xlabel('$x_1$','Interpreter','latex')
ylabel('$\widehat{x}_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Cone Coresponding to X0
patch([X0(1, 1), X0(1, 1), X0(1, 2)], ...
      [X0(1, 2), X0(1, 1), X0(1, 2)], 'r', ...
            'LineWidth', 1, ...
            'FaceAlpha', .7, ...
            'HandleVisibility', 'off');

XE_plot = rect4plot([XE(1:2), XE(3:4)]);
xu_T = XE_plot(:, 1);
xo_T = XE_plot(:, 3);
% Plot Cone Coresponding to Phi^e
patch([xu_T(2, 1), xu_T(2, 1), xo_T(2, 1)], ...
      [xo_T(2, 1), xu_T(2, 1), xo_T(2, 1)], 'g', ...
            'FaceAlpha', .2, ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');

xu = XE_Traj(1:2, :);
xo = XE_Traj(3:4, :);
% Plot trajectory of embedding system that yeilds approximation    
plot(xu(2, :), xo(2, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(2, 1), X0(2, 2),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xu_T(2, 1), xo_T(2, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');

% Plot Diagnol of Embedding Space
plot([-1, 4], [-1, 4], 'k', 'LineWidth', 1.25)

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';


%matlab2tikz('linear_reach5.tikz', 'width', '6cm', 'height', '4cm')

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






