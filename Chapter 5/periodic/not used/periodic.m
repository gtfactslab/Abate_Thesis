%from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000739

clc; clear all;

% Disturbance bound
global W
W = [-.1, .1; ...
     .1, .1];

% Initial Set
X0 = [1, 1.5; ...
      1, 1.5; ...
      1, 1; ...
      0, 0];

% X0 = [-0.0275, 0.1556;
%       -0.1976, 0.0169;
%       0.1, 0.1;
%       -0.7548, -0.7548]

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First part: Simulate Embedding System to Get Donut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_fin = 20;

emb0 = X0(:); % Initial embedding state

[t, xs]   = ode45(@(t,y) F(y, W(:, 2)), [0, T_fin], X0(:, 1));

[te,embs] = ode45(@(t,y) [decomp(y(1:4), W(:, 1), y(5:8), W(:, 2));  ...
                          decomp(y(5:8), W(:, 2), y(1:4), W(:, 1))], ...
                          [0, T_fin],emb0);
         


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on
%axis([-.6, .85, -.4, .5])
axis([-.6, 1.6, -.5, 1.6])
xlabel('$x_1$','Interpreter','latex')
%xticks([-.5,-.25, 0, 0.25, .5, .75])
xticks([-.5, 0, 0.5, 1, 1.5])
ylabel('$x_2$','Interpreter','latex')
%yticks([-.25, 0, .25, .5])
yticks([-.5, 0, .5, 1, 1.5])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Invariant Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('periodic.mat') % Get set
holder = outer(:, 1);
    for i = 2:size(outer, 2)
        if norm(holder(:, end) - outer(:, i)) >= .01
            holder(:, end+1) =  outer(:, i);
        else
        end
    end
patch(holder(1, :), holder(2, :), 'b', 'FaceAlpha', .3, 'LineWidth', 1.3, 'HandleVisibility', 'off');  % Plot

holder = inner(:, 1);
    for i = 2:size(inner, 2)
        if norm(holder(:, end) - inner(:, i)) >= .01
            holder(:, end+1) =  inner(:, i);
        else
        end
    end
patch(inner(1, :), inner(2, :), 'w', 'LineWidth', 1.3, 'HandleVisibility', 'off');                   % Plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get initial points for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = .005;
[a1, a2] = meshgrid(X0(1, 1):dx:X0(1, 2), ...
                    X0(2, 1):dx:X0(2, 2));
a1 = a1(:);
a2 = a2(:);

k = boundary(a1, a2);
a1 = a1(k);
a2 = a2(k);
Initial_Points = zeros(4, size(a1, 1));
for i = 1:1:size(a1, 1)
    Initial_Points(:, i) = [a1(i); a2(i); X0(3, 1); X0(4, 1)];
end
clear a1 a2 k

% plot initial set
patch(Initial_Points(1, :), Initial_Points(2, :), 'r', 'LineWidth', 1.3, 'FaceAlpha', .7)
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduct Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = .005;
T_fin = 14;
dw = .05;

holder = Initial_Points;
holder2 = [];

Embedding_State = X0(:);
corner_holder = [];

for t = 0:dt:T_fin
    t
    
    if t ~= 0
        holder = Next_Set;
    end
    
    if t == .05
        dw = .2;
    end
    
    
    % propegate reachable set forward in time
    for i = 1:size(holder, 2)
        x_now = holder(:, i);
        for w1 = W(1, 1) : dw : W(1, 2)
            for w2 = W(2, 1) : dw : W(2, 2)
                x_next = x_now + dt*F(x_now, [w1, w2]);
                holder2 = [holder2, x_next];
            end
        end
    end
    %k = boundary(holder2(1, :)',holder2(2, :)', fid);
    [k,av] = convhull(holder2([1,2], :)');
    Next_Set = holder2(:, k);
    holder2 = [];

    % propegate embedding state forward in time
    Embedding_State = Embedding_State + dt*E(Embedding_State);

    if t == 10 || t == 12 || t == 14
    if t~=0
        %delete(a)
    end
    rect = makeRectangle([Embedding_State(1:2), Embedding_State(5:6)]);
    patch(rect(1, :), rect(2, :), 'w', 'HandleVisibility', 'off');
    patch(rect(1, :), rect(2, :), 'g', 'FaceAlpha', .2, 'HandleVisibility', 'off');

    holder = Next_Set(:, 1);
    for i = 2:size(Next_Set, 2)
        if norm(holder(:, end) - Next_Set(:, i)) >= .01
            holder(:, end+1) =  Next_Set(:, i);
        else
        end
    end
    patch(holder(1, 1:end), holder(2, 1:end), 'g', 'HandleVisibility', 'off');
    drawnow
    end
end

%%delete(a)
%rect = makeRectangle([Embedding_State(1:2), Embedding_State(5:6)]);
%patch(rect(1, :), rect(2, :), 'm', 'FaceAlpha', .8, 'HandleVisibility', 'off')
%patch(Next_Set(1, 1:2:end), Next_Set(2, 1:2:end), 'g', 'FaceAlpha', .8, 'HandleVisibility', 'off')
%scatter(Embedding_State([1, 5]), Embedding_State([2, 6]), 'k', 'filled', 'HandleVisibility', 'off');
drawnow
Leg = legend();
set(Leg,'visible','off')

%matlab2tikz('periodic.tikz', 'width', '6cm', 'height', '4cm')




function out = F(x, w)
    out = [-2*x(1) + x(2)*(1+x(1)) + x(3) + w(1); ...
           -x(2) + (1 - x(2))*x(1) + w(2); ...
           -x(4); ...
           x(3)];
end

function out = E(a)
    global W
    x = a(1:4);
    xhat = a(5:8);
    
    out = [decomp(x,W(:, 1), xhat, W(:, 2)); ...
           decomp(xhat, W(:, 2), x, W(:, 1))];
end

function out = decomp(x, w, xhat, what)
    if x(1) >= -1
        out(1, 1) = -2*x(1) + x(2)*(1+x(1)) + x(3) + w(1);
    else
        out(1, 1) = -2*x(1) + xhat(2)*(1+x(1)) + x(3) + w(1);
    end
    
    if (1 - x(2)) >= 0 
        out(2, 1) = -x(2) + (1 - x(2))*x(1) + w(2);
    else
        out(2, 1) = -x(2) + (1 - x(2))*xhat(1) + w(2);
    end
    out(3:4, 1) = [-xhat(4); x(3)];
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

    