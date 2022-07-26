clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global A B

A = [0, -1; ...
     1, 0];


X0 = [2, 2.5; ...
      1, 1.5];

B = [1; 0];
W = [0, 1];

T_sim = 4;
dx = .1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: Setup plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; hold on; grid on;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
axis([-3.5, 3, -4, 4])
xticks(-3:3:3);
yticks(-4:4:4);

% plot initial set
X0_plot = rect4plot(X0);
patch(X0_plot(1, :), X0_plot(2, :), 'r', ...
                            'FaceAlpha', .3, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');

scatter(X0(1, :), X0(2, :), 'black', 'filled');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XE, XE_Traj] = PhiE(T_sim, [X0(:, 1); X0(:, 2)], W(:), A, B);


X_initial = rect4plot(X0);




dt = .002;
T_vec = 0:dt:T_sim;
U_sim = W(1):(W(2) - W(1)):W(2);

[X1, X2] = meshgrid(X0(1, 1):(X0(1, 2) - X0(1, 1))/5: X0(1, 2), ... 
                    X0(2, 1):(X0(2, 2) - X0(2, 1))/5: X0(2, 2))

initial_set = [X1(:), X2(:)]';
k = boundary(initial_set(1, :)', initial_set(2, :)');
initial_set = initial_set(:, k);






Reachable_Set = initial_set;
y_now = X0(:);
for i = 1:size(T_vec, 2)-1
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

    t_now = T_vec(i);
    BM = expm(-A*t_now)*B;
    BMp = BM.*(BM >= 0);
    BMm = BM - BMp;
    y_now = y_now + dt*[BMp, BMm; BMm, BMp]*W(:);


    if sum(t_now  == [1, 2, 4] - dt)
        T_now = expm(A*t_now);
        XE_plot = T_now*rect4plot(reshape(y_now, [2, 2]));
        patch(XE_plot(1, :), XE_plot(2, :), [1, 1, 0], ...
                                    'FaceAlpha', .2, ...
                                    'LineWidth', 1.25, ...
                                    'HandleVisibility', 'off');
        points = T_now*[y_now(1:2, 1), y_now(3:4, 1)];
        scatter([points(1, :)], points(2, :), 'black', 'filled')
        
        holder = Reachable_Set(:, 1);
        for i = 2:size(Reachable_Set, 2)
            if norm(holder(:, end) - Reachable_Set(:, i)) >= .02
                holder(:, end+1) =  Reachable_Set(:, i);
            else
            end
        end

        % plot true reachable set
        %patch(holder(1, :), holder(2, :), 'w', 'FaceAlpha', .3)
        patch(holder(1, :), holder(2, :), 'g', ...
                                    'FaceAlpha', .6, ...
                                    'LineWidth', 1.25, ...
                                    'HandleVisibility', 'off')
        drawnow
        
        % plot points

    end
end

% plot MM estimate





Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';



matlab2tikz('linear_reach1.tikz', 'width', '6cm', 'height', '4cm')



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
    dt = .002;
    T_sim = 0:dt:T ;
    Xt = X;
    for i = 1:size(T_sim, 2)-1
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


function dxdt = Dynamics(x, u)
    global A B
    dxdt = A*x + B*u;
end






