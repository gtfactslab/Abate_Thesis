
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 1/13/2021
% Description:  This script generates Figures 1a 1b and 1c.
%               Forward time reachable sets are predicted using MM.

clc; clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Predict Reachable Sets using d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intervals Defining Initial Set
X0 = [-1/2 , 1/2; ....
      -1/2 , 1/2];
% Check to make sure X0 is a valid rectangle
if X0(1, 2) < X0(1, 1) || X0(2, 2) < X0(2, 1)
    print('Error 1')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .005;   % Timestep for Simulation
T  = 1;      % Prediction Time-Horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0_Boundary = makeRectangle(X0);
Phi0 = X0_Boundary;
Phi_size = size(Phi0, 2);

Phi = Phi0;
holder = Phi;
T_size = size(0:dt:T, 2);

xu = zeros(2, T_size + 1);  xu(:, 1) = X0(:, 1);
xo = zeros(2, T_size + 1);  xo(:, 1) = X0(:, 2);
xu3 = zeros(2, T_size + 1);  xu3(:, 1) = X0(:, 1);
xo3 = zeros(2, T_size + 1);  xo3(:, 1) = X0(:, 2);
 
% Compute Time = 1 Second Reachable Set of System
% Compute MM approximation of Reachable Set
W = .2+[0, .1];
load('monotonicity.mat')
x_T = holder;

xo = [.5; .5];
Xo = zeros(2, T_size);
Xo(:, 1) = xo;
xu = -[.5; .5];
Xu = zeros(2, T_size);
Xu(:, 1) = xu;
tic
for i = 1:T_size
    Xu(:, i+1) = Xu(:, i) +dt*dxdt(Xu(:, i), W(:, 1));
    Xo(:, i+1) = Xo(:, i) +dt*dxdt(Xo(:, i), W(:, 2));
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; 
ax = gca;
axis([-1, 5, -1, 3])
xticks([-1, 0, 1, 2, 3, 4, 5])
yticks([-1, 0, 1, 2, 3])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')


xu_T = Xu(:, end);
xo_T = Xo(:, end);

rect = [xu_T(1, 1), xo_T(1, 1), xo_T(1, 1), xu_T(1, 1); ...
        xu_T(2, 1), xu_T(2, 1), xo_T(2, 1), xo_T(2, 1)]
% Plot Initial Set X0
patch(rect(1, :), rect(2, :), 'r' , ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .1, ...
            'HandleVisibility', 'off');
        
% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r' , ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .1, ...
            'HandleVisibility', 'off');
% Plot MM Overapproximation of RF(1, X0) with d2
                  

% Plot Time = 1 Reachable Set RF(1, X0)             
k = boundary(x_T(1, :)',x_T(2, :)',.8);
patch(x_T(1, k), x_T(2, k), 'g', ...
                            'FaceAlpha', .1, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');

plot(Xu(1, :), Xu(2, :), 'b')
plot(Xo(1, :), Xo(2, :), 'b')
                        
Leg = legend();
set(Leg,'visible','off')

%matlab2tikz('monotonicity.tikz', 'width', '4cm', 'height', '3cm')












%%






figure(2); clf;
hold on; 
ax = gca;
axis([-1.5, 1.5, -1.1, 0])
xticks([-1, 0, 1, 2, 3, 4, 5])
yticks([-1, 0, 1, 2, 3])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')


xo = -1+[.4; .4];
Xo = zeros(2, T_size);
Xo(:, 1) = xo;
xu = -1+[0.1; .3];
Xu = zeros(2, T_size);
Xu(:, 1) = xu;
tic
for i = 1:T_size
    Xu(:, i+1) = Xu(:, i) +dt*dxdt(Xu(:, i), W(:, 1));
    Xo(:, i+1) = Xo(:, i) +dt*dxdt(Xo(:, i), W(:, 2));
end

plot(Xu(1, :), Xu(2, :), 'b')
plot(Xo(1, :), Xo(2, :), 'b')


for i = 1:60:T_size
   scatter([Xu(1, i) Xo(1, i)], [Xu(2, i) Xo(2, i)], 'filled');
    plot([Xu(1, i), Xu(1, i), Xu(1, i)+.11], ...
         [Xu(2, i)+.06, Xu(2, i), Xu(2, i)], 'r')
end

                        
Leg = legend();
set(Leg,'visible','off')

matlab2tikz('monotonicity5.tikz', 'width', '4cm', 'height', '3cm')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/20 :X0(1, 2), ...
                            X0(2, 1): d(2)/20 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

function out = dxdt(x, w)
    out = [ x(2, 1)^3 - x(1, 1)+2+w;...
            x(1, 1)];
end

