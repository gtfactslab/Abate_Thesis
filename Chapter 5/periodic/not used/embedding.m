%from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000739

clc; clear all;

% Disturbance bound
global W
W = [-.1, .1; ...
     .1, .1];

% Initial Set
X0 = [1, 1.5; ...
      1, 1.5; ...
      1, 1.0; ...
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
         


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on

subplot(4, 1, 1); hold on;
%axis([-.6, .85, -.4, .5])
%axis([-.6, 1.6, -.5, 1.6])
%xticks([-.5,-.25, 0, 0.25, .5, .75])
%xticks([-.5, 0, 0.5, 1, 1.5])
%yticks([-.25, 0, .25, .5])
%yticks([-.5, 0, .5, 1, 1.5])
xlabel('$t$','Interpreter','latex')
ylabel('$x_1$','Interpreter','latex')
plot(te', embs(:, 1)', 'b')
plot(te', embs(:, 5)', 'r')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

subplot(4, 1, 2); hold on;
xlabel('$t$','Interpreter','latex')
ylabel('$x_1$','Interpreter','latex')
plot(te', embs(:, 2)', 'b')
plot(te', embs(:, 6)', 'r')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

subplot(4, 1, 3); hold on;
xlabel('$t$','Interpreter','latex')
ylabel('$x_1$','Interpreter','latex')
plot(te', embs(:, 3)', 'b')
plot(te', embs(:, 7)', 'r')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

subplot(4, 1, 4); hold on;
xlabel('$t$','Interpreter','latex')
ylabel('$x_1$','Interpreter','latex')
plot(te', embs(:, 4)', 'b')
plot(te', embs(:, 8)', 'r')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

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

    