
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 1/13/2021
% Description:  This script generates Figure 2.
%               An invariant rectangle is computed via analysis of 
%               equilibria in the embedding space.  The smallest robustly 
%               forward invariant set in X is also computed via exhaustive
%               simulation.
%               This script requires the external file Ex2_data.mat,
%               availible though the GaTech FactsLab GitHub

clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global W
W = [-1, 1]; % Didurbance Bound


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Hyperrectangular RFI Set Using MM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute an equilibrium in embedding space.
% This point defines a robustly forward invairant set for the original
% dynamics.
tic
thing = fsolve(@E, 1*[-1; -1; 1; 1]);
toc
XE = reshape(thing, 2, 2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Smallest Attractive RFI Set Via Exhaustive Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an exhaustive simulation was computed externally, and the resulting 
% reachable set is provided in the external file Ex2_data.mat
XE_Boundary = makeRectangle(XE);
load('Ex2_data.mat')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-2, 1.5, -1.5, 4]);
xticks([-2 -1 0 1 ])
yticks([-1 0 1 2 3 4])

% Plot RFI set from MM approach
patch(XE_Boundary(1, :), XE_Boundary(2, :), 'r', ...
            'LineWidth', 1.15, ...
            'FaceAlpha', .2, ...
            'HandleVisibility', 'off');
scatter(XE(1, :), XE(2, :), 'k', 'filled', ...
            'HandleVisibility', 'off');
% Plot RFI set from exhaustive simulation
patch(REACH(1, 1:8:end), REACH(2, 1:8:end), 'g', ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .8, ...
            'HandleVisibility', 'off');

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

%matlab2tikz('invariance_1.tikz', 'width', '3.5cm', 'height', '2cm')






%%

figure(2); clf;
hold on; grid on;
xlabel('$x_1$','Interpreter','latex')
ylabel('$\widehat{x}_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-3, 3, -3, 3]);
xticks([-2 -1 0 1, 2 ])
yticks([-2 -1 0 1 2])

plot(e_holder1(1, :), e_holder1(3, :))
scatter(e_holder1(1, end), e_holder1(3, end))

plot([-3, 3], [-3, 3])
Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

e0 = [-2.5; -1; -1; 0];
e_holder1 = e0;
for i = 1:2000
    e_holder1(:, i+1) = e_holder1(:, i) + .01*E(e_holder1(:, i));
end
plot(e_holder1(1, :), e_holder1(3, :))


e0 = [1; -2; 2; -1];
e_holder1 = e0;
for i = 1:2000
    e_holder1(:, i+1) = e_holder1(:, i) + .01*E(e_holder1(:, i));
end
plot(e_holder1(1, :), e_holder1(3, :))


%matlab2tikz('invariance_2.tikz', 'width', '3.5cm', 'height', '2cm')



%%

figure(3); clf;
hold on; grid on;
xlabel('$x_1$','Interpreter','latex')
ylabel('$\widehat{x}_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-3, 0, -3, 1]);
%xticks([-2 -1 0 1, 2 ])
%yticks([-2 -1 0 1 2])

plot([-3, 3], [-3, 3], 'k',  'LineWidth', 1.2)
Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

q = 70;
e0 = [-2.25; -1; -.5; -1];
e_holder1 = e0;

e0 = [-2.75; -1; 0; -1];
e_holder2 = e0;
for i = 1:349
    e_holder1(:, i+1) = e_holder1(:, i) + .001*E(e_holder1(:, i));
    e_holder2(:, i+1) = e_holder2(:, i) + .001*E(e_holder2(:, i));
end
plot(e_holder1(1, :), e_holder1(3, :), 'b')

plot(e_holder2(1, :), e_holder2(3, :), 'r')

for i = 1:70:349
    thing1 = [e_holder2(1, i), e_holder2(1, i), e_holder2(1, i) + .17];
    thing2 = [e_holder2(3, i)-.25, e_holder2(3, i), e_holder2(3, i)];
    plot(thing1, thing2, 'r')
    
    scatter([e_holder1(1, i), e_holder2(1, i)], [e_holder1(3, i), e_holder2(3, i)], 'filled')
end

matlab2tikz('emb_mon.tikz', 'width', '3.5cm', 'height', '2cm')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = F(x, w)
    out = [ - x(1)^3 - x(2) - w; ...
            x(1)^2 - x(2) + w^3]; 
end


function out = d(x, w, xhat, what)
    out(1, 1) = -x(1)^3 -xhat(2) - what;
    
    if x(1, 1) >= 0 && x(1, 1) >= - xhat(1, 1)  
        out(2, 1) = x(1)^2 - x(2) + w^3;   
    elseif xhat(1, 1) <= 0 && x(1, 1) <= -xhat(1, 1)
        out(2, 1) = xhat(1)^2 - x(2) + w^3;  
    elseif (x(1, 1) <= 0) && (xhat(1, 1) >= 0)
        out(2, 1) = - x(2) + w^3;
    end
end

function out = E(in)
    global W
    x = in(1:2);
    xhat = in(3:4);
    
    out = [d(x, W(:, 1), xhat, W(:, 2)); ...
           d(xhat, W(:, 2), x, W(:, 1))];
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end
