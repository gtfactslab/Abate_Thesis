%from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000739

clc; clear all;

% Disturbance bound
global W
W = [-.1, .1; ...
     .1, .1];

% Initial Set
X0 = [1, 1.5; ...
      1, 1.5; ...
      1, 1.5; ...
      0, 0.5];

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First part: Simulate Embedding System to Get Donut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_fin = 30;

emb0 = X0(:); % Initial embedding state

opts = odeset('RelTol',1e-4,'AbsTol',1e-5);
[t, xs]   = ode45(@(t,y) F(y, W(:, 2)), [0, T_fin], X0(:, 1));

[te,embs] = ode45(@(t,y) [decomp(t, y(1:4), W(:, 1), y(5:8), W(:, 2));  ...
                          decomp(t, y(5:8), W(:, 2), y(1:4), W(:, 1))], ...
                          [0, T_fin],emb0, opts);
         


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on

embs = embs';
te = te';
embs = embs(:, 1:5:end);
te = te(:, 1:5:end);

subplot(4, 1, 1); hold on; grid on
%xlabel('$t$','Interpreter','latex')
ylabel('$y_1$','Interpreter','latex')
axis([0, T_fin, -1, 2])
yticks(-1:2);
patch([te,  fliplr(te)], [embs(1, :), fliplr(embs(5, :))], 'b', 'FaceAlpha', .1)
plot(te, embs(1, :), 'r', 'LineWidth', 1.5)
plot(te, embs(5, :), 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
drawnow
Leg = legend();
set(Leg,'visible','off')


subplot(4, 1, 2); hold on; grid on
axis([0, T_fin, -1.1, 2])
yticks(-1:2);
%xlabel('$t$','Interpreter','latex')
ylabel('$y_2$','Interpreter','latex')
patch([te,  fliplr(te)], [embs(2, :), fliplr(embs(6, :))], 'b', 'FaceAlpha', .1)
plot(te, embs(2, :), 'r', 'LineWidth', 1.5)
plot(te, embs(6, :), 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
drawnow
Leg = legend();
set(Leg,'visible','off')


subplot(4, 1, 3); hold on; grid on
axis([0, T_fin, 0.5, 2])
yticks(1:2);
%xlabel('$t$','Interpreter','latex')
ylabel('$y_3$','Interpreter','latex')
patch([te,  fliplr(te)], [embs(3, :), fliplr(embs(7, :))], 'b', 'FaceAlpha', .1)
plot(te, embs(3, :), 'r', 'LineWidth', 1.5)
plot(te, embs(7, :), 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
drawnow
Leg = legend();
set(Leg,'visible','off')


subplot(4, 1, 4); hold on; grid on
axis([0, T_fin, -0.5, 1])
yticks(0:1);
xlabel('$t$','Interpreter','latex')
ylabel('$y_4$','Interpreter','latex')
patch([0, T_fin,  T_fin, 0], [embs(4, 1), embs(4, end), embs(8, end), embs(8, 1)], 'b', 'FaceAlpha', .1)
plot([0, T_fin], [embs(4, 1), embs(4, end)], 'r', 'LineWidth', 1.5)
plot([0, T_fin], [embs(8, 1), embs(8, end)], 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
%
drawnow
Leg = legend();
set(Leg,'visible','off')

% matlab2tikz('periodic.tikz', 'width', '6cm', 'height', '4cm')

holder = [];
for i = 1:size(te, 2)
    t = te(i);
    T = Transform(t);
    Tp = (T >= 0).*T;
    Tm = (T <= 0).*T;
    TT = [Tp, Tm; Tm, Tp];

    state1 = inv(TT)*embs(:, i);

    holder = [holder, [state1]];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); clf;
hold on; grid on

embs = holder;
subplot(4, 1, 1); hold on; grid on
%xlabel('$t$','Interpreter','latex')
ylabel('$y_1$','Interpreter','latex')
axis([0, T_fin, -1, 2])
yticks(-1:2);
patch([te,  fliplr(te)], [embs(1, :), fliplr(embs(5, :))], 'b', 'FaceAlpha', .1)
plot(te, embs(1, :), 'r', 'LineWidth', 1.5)
plot(te, embs(5, :), 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
drawnow
Leg = legend();
set(Leg,'visible','off')


subplot(4, 1, 2); hold on; grid on
axis([0, T_fin, -1.1, 2])
yticks(-1:2);
%xlabel('$t$','Interpreter','latex')
ylabel('$y_2$','Interpreter','latex')
patch([te,  fliplr(te)], [embs(2, :), fliplr(embs(6, :))], 'b', 'FaceAlpha', .1)
plot(te, embs(2, :), 'r', 'LineWidth', 1.5)
plot(te, embs(6, :), 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
drawnow
Leg = legend();
set(Leg,'visible','off')


subplot(4, 1, 3); hold on; grid on
axis([0, T_fin, 0.5, 2])
yticks(1:2);
%xlabel('$t$','Interpreter','latex')
ylabel('$y_3$','Interpreter','latex')
patch([te,  fliplr(te)], [embs(3, :), fliplr(embs(7, :))], 'b', 'FaceAlpha', .1)
plot(te, embs(3, :), 'r', 'LineWidth', 1.5)
plot(te, embs(7, :), 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
drawnow
Leg = legend();
set(Leg,'visible','off')


subplot(4, 1, 4); hold on; grid on
axis([0, T_fin, -0.5, 1])
yticks(0:1);
xlabel('$t$','Interpreter','latex')
ylabel('$y_4$','Interpreter','latex')
patch([0, T_fin,  T_fin, 0], [embs(4, 1), embs(4, end), embs(8, end), embs(8, 1)], 'b', 'FaceAlpha', .1)
plot([0, T_fin], [embs(4, 1), embs(4, end)], 'r', 'LineWidth', 1.5)
plot([0, T_fin], [embs(8, 1), embs(8, end)], 'b', 'LineWidth', 1.5)
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
%
drawnow
Leg = legend();
set(Leg,'visible','off')




function out = Transform(t)
    out = [1, 0, 0, 0; ...
           0, 1, 0, 0; ...
           0, 0, cos(t), -sin(t); ...
           0, 0, sin(t),  cos(t)];
end

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

function out = decomp(t, x, w, xhat, what)
    thing = max([cos(t), 0])*x(3) + min([cos(t), 0])*xhat(3) + ...
            max([sin(t), 0])*x(4) + min([sin(t), 0])*xhat(4);
    if x(1) >= -1
        out(1, 1) = -2*x(1) + x(2)*(1+x(1)) + thing + w(1);
    else
        out(1, 1) = -2*x(1) + xhat(2)*(1+x(1)) + thing + w(1);
    end
    
    if (1 - x(2)) >= 0 
        out(2, 1) = -x(2) + (1 - x(2))*x(1) + w(2);
    else
        out(2, 1) = -x(2) + (1 - x(2))*xhat(1) + w(2);
    end
    out(3:4, 1) = [0; 0];
end



function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

    