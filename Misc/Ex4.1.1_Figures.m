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
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First part: Simulate Embedding System to Get Donut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_fin = 20;

emb0 = X0(:); % Initial embedding state

[t, xs]   = ode45(@(t,y) F(y, W(:, 2)), [0, T_fin], X0(:, 1));

[te,embs] = ode45(@(t,y) [decomp(y(1:4), W(:, 1), y(5:8), W(:, 2));  ...
                          decomp(y(5:8), W(:, 2), y(1:4), W(:, 1))], ...
                          [0, T_fin],emb0);
         
% Figure 1
figure(1); clf;

% Subplot 1
ax=subplot(2,1,1);
hold on; grid on;
xlabel('$t$','Interpreter','latex')
ylabel('$x_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([0, T_fin, -.6, 2])

fill([te;flipud(te)],[embs(:,1);flipud(embs(:,5))],'g','FaceAlpha',.1, 'HandleVisibility', 'off');
plot(t,xs(:,1), 'b', 'LineWidth', 2, 'HandleVisibility', 'off');

Leg = legend();
set(Leg,'visible','off')

% Subplot 2
ax=subplot(2,1,2);
hold on; grid on;
xlabel('$t$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([0, T_fin, -.5, 1.5])

fill([te;flipud(te)],[embs(:,2);flipud(embs(:,6))],'g','FaceAlpha',.1, 'HandleVisibility', 'off');
plot(t,xs(:,2), 'b', 'LineWidth', 2, 'HandleVisibility', 'off');

Leg = legend();
set(Leg,'visible','off')

matlab2tikz('F3a.tikz', 'width', '6cm', 'height', '4cm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); clf;
hold on; grid on
axis([-.6, 1.6, -.5, 1.6])
xlabel('$x_1$','Interpreter','latex')
xticks([-.5, 0, .5, 1, 1.5])
ylabel('$x_2$','Interpreter','latex')
yticks([-.5, 0, .5, 1, 1.5])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Invariant Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Ex4_Data.mat') % Get set
patch(outer(1, :), outer(2, :), 'b', 'FaceAlpha', .2, 'HandleVisibility', 'off');  % Plot
patch(inner(1, :), inner(2, :), 'w', 'HandleVisibility', 'off');                   % Plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get initial points for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = .001;
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
patch(Initial_Points(1, :), Initial_Points(2, :), 'r')
scatter(X0(1, :), X0(2, :), 'k', 'filled', 'HandleVisibility', 'off');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduct Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = .002;
T_fin = 12;
dw = .05;
fid = .7;

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
        dw = .1;
        fid = 0;
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
    k = boundary(holder2(1, :)',holder2(2, :)', fid);
    Next_Set = holder2(:, k);
    holder2 = [];

    % propegate embedding state forward in time
    Embedding_State = Embedding_State + dt*E(Embedding_State);
    if t~=0
        delete(a)
    end
    rect = makeRectangle([Embedding_State(1:2), Embedding_State(5:6)]);
    a(1) = patch(rect(1, :), rect(2, :), 'm', 'FaceAlpha', .8, 'HandleVisibility', 'off');
    a(2) = patch(Next_Set(1, 1:2:end), Next_Set(2, 1:2:end), 'g', 'HandleVisibility', 'off');
    a(3) = scatter(Embedding_State([1, 5]), Embedding_State([2, 6]), 'k', 'filled', 'HandleVisibility', 'off');
    drawnow
end

delete(a)
rect = makeRectangle([Embedding_State(1:2), Embedding_State(5:6)]);
patch(rect(1, :), rect(2, :), 'm', 'FaceAlpha', .8, 'HandleVisibility', 'off')
patch(Next_Set(1, 1:2:end), Next_Set(2, 1:2:end), 'g', 'FaceAlpha', .8, 'HandleVisibility', 'off')
scatter(Embedding_State([1, 5]), Embedding_State([2, 6]), 'k', 'filled', 'HandleVisibility', 'off');

Leg = legend();
set(Leg,'visible','off')

matlab2tikz('F3b.tikz', 'width', '6cm', 'height', '4cm')




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

    