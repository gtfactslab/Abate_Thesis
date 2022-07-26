
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
dt = .002;   % Timestep for Simulation
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

xu2 = zeros(2, T_size + 1);  xu2(:, 1) = X0(:, 1);
xo2 = zeros(2, T_size + 1);  xo2(:, 1) = X0(:, 2);
xu3 = zeros(2, T_size + 1);  xu3(:, 1) = X0(:, 1);
xo3 = zeros(2, T_size + 1);  xo3(:, 1) = X0(:, 2);
 
% Compute Time = 1 Second Reachable Set of System
% Compute MM approximation of Reachable Set
for t = 1:T_size
    % reachable set computation
    holder2 = zeros(2, Phi_size);
    for i = 1:size(holder, 2)
            x = holder(:, i);
            x_next = x + dt*dxdt(x);
            holder2(:, i) = x_next;
    end
    holder = holder2;
    Phi = [Phi, holder2];
    
    % MM approximation
    xu2(:, t + 1) = xu2(:, t) + dt*d2(xu2(:, t), xo2(:, t));
    xo2(:, t + 1) = xo2(:, t) + dt*d2(xo2(:, t), xu2(:, t));
    xu3(:, t + 1) = xu3(:, t) + dt*d3(xu3(:, t), xo3(:, t));
    xo3(:, t + 1) = xo3(:, t) + dt*d3(xo3(:, t), xu3(:, t));
end

x_T = holder;

xu_T2 = xu2(:, t + 1);
xo_T2 = xo2(:, t + 1);

xu_T3 = xu3(:, t + 1);
xo_T3 = xo3(:, t + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
ax = gca;
axis([-1, 5, -1, 3])
xticks([-1, 0, 1, 2, 3, 4, 5])
yticks([-1, 0, 1, 2, 3])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r' , ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .7, ...
            'HandleVisibility', 'off');
% Plot MM Overapproximation of RF(1, X0) with d2
SF = makeRectangle([xu_T2, xo_T2]);
patch(SF(1, :), SF(2, :), 'm', ...
                          'FaceAlpha', .1, ...
                          'LineWidth', 1.25, ...
                          'HandleVisibility', 'off');
% Plot MM Overapproximation of RF(1, X0) with d3
SF = makeRectangle([xu_T3, xo_T3]);
patch(SF(1, :), SF(2, :), 'g', ...
                          'FaceAlpha', .1, ...
                          'LineWidth', 1.25, ...
                          'HandleVisibility', 'off');
                      

% Plot Time = 1 Reachable Set RF(1, X0)
patch(x_T(1, :), x_T(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');
scatter(xu_T2(1, 1), xu_T2(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(xo_T2(1, 1), xo_T2(2, 1), 'k', 'filled', 'HandleVisibility', 'off');

scatter(xu_T3(1, 1), xu_T3(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(xo_T3(1, 1), xo_T3(2, 1), 'k', 'filled', 'HandleVisibility', 'off');


Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';
%matlab2tikz('F1b.tikz', 'width', '6cm', 'height', '4cm')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;
hold on; grid on;
ax = gca;
axis([-1, 5, -1, 5])
xticks([-1, 0, 1, 2, 3, 4, 5])
yticks([-1, 0, 1, 2, 3, 4, 5])
xlabel('$x_1$','Interpreter','latex')
ylabel('$\widehat{x}_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Cone Coresponding to X0
patch([X0(1, 1), X0(1, 1), X0(1, 2)], ...
      [X0(1, 2), X0(1, 1), X0(1, 2)], 'r', ...
            'LineWidth', 1, ...
            'FaceAlpha', .7, ...
            'HandleVisibility', 'off');
% Plot Cone Coresponding to Phi^e2
patch([xu_T2(1, 1), xu_T2(1, 1), xo_T2(1, 1)], ...
      [xo_T2(1, 1), xu_T2(1, 1), xo_T2(1, 1)], 'm', ...
            'FaceAlpha', .2, ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');
% Plot trajectory of embedding system that yeilds approximation    
plot(xu2(1, :), xo2(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 1), X0(1, 2),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xu_T2(1, 1), xo_T2(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');

            
% Plot Cone Coresponding to Phi^e3
patch([xu_T3(1, 1), xu_T3(1, 1), xo_T3(1, 1)], ...
      [xo_T3(1, 1), xu_T3(1, 1), xo_T3(1, 1)], 'g', ...
            'FaceAlpha', .2, ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');
% Plot trajectory of embedding system that yeilds approximation    
plot(xu3(1, :), xo3(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 1), X0(1, 2),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xu_T3(1, 1), xo_T3(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');

% Plot Diagnol of Embedding Space
plot([-1, 5], [-1, 5], 'k', 'LineWidth', 1.25)

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%matlab2tikz('F1c.tikz', 'width', '6cm', 'height', '4cm')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Compare d and d'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = .4; % Prediction Time-Horizon

X0_Boundary = makeRectangle(X0);
Phi0 = X0_Boundary;
Phi_size = size(Phi0, 2);

Phi = Phi0;
holder = Phi;
T_size = size(0:dt:T, 2);
xu1(:, 1) = X0(:, 1);
xo1(:, 1) = X0(:, 2);
xu2(:, 1) = X0(:, 1);
xo2(:, 1) = X0(:, 2);
xu3(:, 1) = X0(:, 1);
xo3(:, 1) = X0(:, 2);

% Compute Time = .25 Second Reachable Set of System
% Compute two MM approximations of Reachable Set
flag = 0;
t = 0;
while flag == 0
    t = t + 1;
    % get next reachable set
    holder2 = zeros(2, Phi_size);
    for i = 1:size(holder, 2)
            x = holder(:, i);
            x_next = x + dt*dxdt(x);
            holder2(:, i) = x_next;
    end
    holder = holder2;
    Phi = [Phi, holder2];
    
    % MM approximation 1: decompostion function from special case
    xu1(:, t + 1) = xu1(:, t) + dt*d1(xu1(:, t), xo1(:, t));
    xo1(:, t + 1) = xo1(:, t) + dt*d1(xo1(:, t), xu1(:, t));
    
    % MM approximation 2: decompostion function from polynomials
    xu2(:, t + 1) = xu2(:, t) + dt*d2(xu2(:, t), xo2(:, t));
    xo2(:, t + 1) = xo2(:, t) + dt*d2(xo2(:, t), xu2(:, t));
    
    % MM approximation 1: tight decompostion function
    xu3(:, t + 1) = xu3(:, t) + dt*d3(xu3(:, t), xo3(:, t));
    xo3(:, t + 1) = xo3(:, t) + dt*d3(xo3(:, t), xu3(:, t));
    
    if ~(prod( xu1(:, t + 1) >= -5 ) &&  prod( xo1(:, t + 1) <= 5 )  )
        flag = 1;
    end
end
x_T = holder;
xu_T1 = xu1(:, t + 1);
xo_T1 = xo1(:, t + 1);
xu_T2 = xu2(:, t + 1);
xo_T2 = xo2(:, t + 1);
xu_T3 = xu3(:, t + 1);
xo_T3 = xo3(:, t + 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf;
hold on; grid on;
ax = gca;
axis([-4, 6, -1.25, 1.25])
xticks([-4, -2, 0, 2, 4, 6])
yticks([-1, 0, 1])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Prediction From d1
SF1 = makeRectangle([xu_T1, xo_T1]);
patch(SF1(1, :), SF1(2, :), 'b', ...
                            'FaceAlpha', .2, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');
% Create whitespace
SF3 = makeRectangle([xu_T3, xo_T3]);
% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r', ...
                          'FaceAlpha', .7, ...
                          'LineWidth', 1.25, ...
                          'HandleVisibility', 'off');      

                      
% Plot Prediction from d2
SF2 = makeRectangle([xu_T2, xo_T2]);
patch(SF2(1, :), SF2(2, :), 'm', ...
                          'FaceAlpha', .2, ...
                          'LineWidth', 1.25, ...
                          'HandleVisibility', 'off');
                      
% Plot Prediction from d3
%SF3 = makeRectangle([xu_T3, xo_T3]);
%patch(SF3(1, :), SF3(2, :), 'g', ...
%                          'FaceAlpha', .2, ...
%                          'LineWidth', 1.25, ...
%                          'HandleVisibility', 'off');
                      
% Plot Reachable set RF(1/2, X0)
patch(x_T(1, :), x_T(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');
                        

% Dot the relevant rectangle corners
scatter(xu_T1(1, 1), xu_T1(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(xo_T1(1, 1), xo_T1(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(xu_T2(1, 1), xu_T2(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(xo_T2(1, 1), xo_T2(2, 1), 'k', 'filled', 'HandleVisibility', 'off');

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%matlab2tikz('F1a.tikz', 'width', '6cm', 'height', '4cm')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = d1(x, xh)
    out = [xh(2)^2 + 2 + 10*(x(2) - xh(2)); ...
           x(1)]; 
end

function out = d2(x, xh)
    if x(2, 1) >= 0 && x(2, 1) >= - xh(2, 1)  
        out(1, 1) = x(2, 1)^2 + 2;   
    elseif xh(2, 1) <= 0 && x(2, 1) <= -xh(2, 1)
        out(1, 1) = xh(2, 1)^2 +  2;  
    elseif (x(2, 1) <= 0) && (xh(2, 1) >= 0)
        out(1, 1) = x(2, 1)*xh(2, 1) + 2;
    end

    out(2, 1) = x(1); 
end

function out = d3(x, xh)
    if x(2, 1) >= 0 && x(2, 1) >= - xh(2, 1)  
        out(1, 1) = x(2, 1)^2 + 2;   
    elseif xh(2, 1) <= 0 && x(2, 1) <= -xh(2, 1)
        out(1, 1) = xh(2, 1)^2 +  2;  
    elseif (x(2, 1) <= 0) && (xh(2, 1) >= 0)
        out(1, 1) = 2;
    end

    out(2, 1) = x(1); 
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/10 :X0(1, 2), ...
                            X0(2, 1): d(2)/10 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

function out = dxdt(x)
    out = [ x(2, 1)^2 + 2;...
            x(1, 1) ];
end