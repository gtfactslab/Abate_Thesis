clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global A B

V = [ 2,  1; ...
       1, 2];

L = [-1.5, 0; ...
     0, .1];

A = V*L*inv(V);


B = [3; 1];

AM = [A, B; zeros(1, 3)];
X0 = [1; ...
      0];

W = [1, 2];

T_sim = 1;
dx = .1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XE, XE_Traj] = PhiE(T_sim, [X0; X0], W(:), A, B);


X_initial = X0 ;

Reachable_Set = Reach(T_sim, X_initial, W);
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
%axis([-0.5, 1.75, -0.5, 4.5])
%xticks(-1:4);
%yticks(-1:3);

% plot initial set
scatter(X0(1, 1), X0(2, 1), 80, 'r', 'filled')
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



[XE, XE_Traj] = PhiE(T_sim, [inv(V)*X0; inv(V)*X0], W(:), L, inv(V)*B);
XE_plot = V*rect4plot([XE(1), XE(3); XE(2), XE(4)]);
patch(XE_plot(1, :), XE_plot(2, :), 'r', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');

%%

[V, L] = eig(AM);
[XE, XE_Traj] = PhiE(T_sim, [inv(V)*[X0; W(1)]; inv(V)*[X0; W(2)]], [0; 0], L, [0;0;0]);
%XE_plot = V(1:2, 1:2)*rect4plot([XE(1), XE(4); XE(2), XE(5)]);

clc
XE = reshape(XE, [3, 2])

holder = [];
for i = 1:2
    for j = 1:2
        for k = 1:2
            holder = [holder, V*[XE(1, i); XE(2, j);XE(3, k)] ];
        end
    end
end


scatter(holder(1, :), holder(2, :), 'b')
% patch(XE_plot(1, :), XE_plot(2, :), 'g', ...
%                             'FaceAlpha', .3, ...
%                             'LineWidth', 1.25, ...
%                             'HandleVisibility', 'off');


%matlab2tikz('linear_reach5.tikz', 'width', '6cm', 'height', '4cm')

function out = rect4plot(X)
    out = [X(1, 1), X(1, 2), X(1, 2), X(1, 1); ...
           X(2, 1), X(2, 1), X(2, 2), X(2, 2)];
end

function out = Decomposition(X, U, A, B)
    Ap =  (A - diag(diag(A))).*(A >= 0) + diag(diag(A));
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
    Xflip = [X(end/2+1:end, 1); X(1:end/2, 1)];
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
        if i >= 2
            k = convhull(X_Next(1, :)', X_Next(2, :)' );
            Reachable_Set = X_Next(:, k);
        else
            Reachable_Set = X_Next;
        end
    end
end

function dxdt = Dynamics(x, u)
    global A B
    dxdt = A*x + B*u;
end






