
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 9/9/2021
% Description:  This script generates Figures 1a and 1b.
%               Forward time reachable sets are predicted using MM.

clc; clear all; %close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Predict Reachable Sets using d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intervals Defining Initial Set
X0 = [-1/2 , 1/2];
W = [-pi/8, pi/4];

X0 = [X0; W];

% Check to make sure X0 is a valid rectangle
if X0(1, 2) < X0(1, 1) || X0(2, 2) < X0(2, 1)
    print('Error 1')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .002;   % Timestep for Simulation
T  = 2;      % Prediction Time-Horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0_Boundary = makeRectangle(X0);
Phi0 = X0_Boundary;
Phi_size = size(Phi0, 2);

Phi = Phi0;
holder = Phi;
T_size = size(0:dt:T, 2);
xu1 = zeros(2, T_size + 1);  xu1(:, 1) = X0(:, 1);
xo1 = zeros(2, T_size + 1);  xo1(:, 1) = X0(:, 2);
xu2 = zeros(2, T_size + 1);  xu2(:, 1) = X0(:, 1);
xo2 = zeros(2, T_size + 1);  xo2(:, 1) = X0(:, 2);
xu3 = zeros(2, T_size + 1);  xu3(:, 1) = X0(:, 1);
xo3 = zeros(2, T_size + 1);  xo3(:, 1) = X0(:, 2);
xu4 = zeros(2, T_size + 1);  xu4(:, 1) = X0(:, 1);
xo4 = zeros(2, T_size + 1);  xo4(:, 1) = X0(:, 2);
xu5 = zeros(2, T_size + 1);  xu5(:, 1) = X0(:, 1);
xo5 = zeros(2, T_size + 1);  xo5(:, 1) = X0(:, 2);


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
    xu1(:, t + 1) = xu1(:, t) + dt*d1(xu1(:, t), xo1(:, t));
    xo1(:, t + 1) = xo1(:, t) + dt*d1(xo1(:, t), xu1(:, t));

    % MM approximation 2
    xu2(:, t + 1) = xu2(:, t) + dt*d2(xu2(:, t), xo2(:, t));
    xo2(:, t + 1) = xo2(:, t) + dt*d2(xo2(:, t), xu2(:, t));

     % MM approximation 2
    xu3(:, t + 1) = xu3(:, t) + dt*d3(xu3(:, t), xo3(:, t));
    xo3(:, t + 1) = xo3(:, t) + dt*d3(xo3(:, t), xu3(:, t));

    % MM approximation 2
    xu4(:, t + 1) = xu4(:, t) + dt*d4(xu4(:, t), xo4(:, t));
    xo4(:, t + 1) = xo4(:, t) + dt*d4(xo4(:, t), xu4(:, t));

    % MM approximation 2
    xu5(:, t + 1) = xu5(:, t) + dt*d5(xu5(:, t), xo5(:, t));
    xo5(:, t + 1) = xo5(:, t) + dt*d5(xo5(:, t), xu5(:, t));
end
x_T = holder;
xu_T1(:, 1) = xu1(:, t + 1);
xu_T1(:, 2) = xu2(:, t + 1);
xu_T1(:, 3) = xu3(:, t + 1);
xu_T1(:, 4) = xu4(:, t + 1);
xu_T1(:, 5) = xu5(:, t + 1);

xo_T1(:, 1) = xo1(:, t + 1);
xo_T1(:, 2) = xo2(:, t + 1);
xo_T1(:, 3) = xo3(:, t + 1);
xo_T1(:, 4) = xo4(:, t + 1);
xo_T1(:, 5) = xo5(:, t + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
ax = gca;
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-3*pi/4, 7*pi/4, -1, 3])
xticks([-pi, -pi/2, 0, pi/2, pi, 3*pi/2])
yticks([-1:6])
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$','$3\pi/2$'})
yticklabels({'','$X$','$d^{1}(w, \widehat{w})$', '$d^{2}(w, \widehat{w})$'})
xlabel('$x$','Interpreter','latex')


% Plot Initial Set X0
plot([min(X0_Boundary(1, :)), max(X0_Boundary(1, :))], [0, 0], 'r' , ...
            'LineWidth', 2, ... 
            'HandleVisibility', 'off');
plot([min(X0_Boundary(1, :)), min(X0_Boundary(1, :))], [-.1, +.1], 'r' , ...
                'LineWidth', 2, ... 
                'HandleVisibility', 'off');
plot([max(X0_Boundary(1, :)), max(X0_Boundary(1, :))], [-.1, +.1], 'r' , ...
                'LineWidth', 2, ... 
                'HandleVisibility', 'off');

c = [0, .8, 0];
% Plot Initial Set reachable set 1 Plot Time = 1 Reachable Set RF(1, X0)
plot([min(x_T(1, :)), max(x_T(1, :))], [0, 0], 'Color', c, ...
                            'LineWidth', 2, ...
                            'HandleVisibility', 'off');

% Plot Initial Set reachable set 1 Plot Time = 1 Reachable Set RF(1, X0)
plot([min(x_T(1, :)), min(x_T(1, :))], [-1, 6], 'k--', ...
                            'Color', [0.6, 0.6, 0.6], ...
                            'LineWidth', 2, ...
                            'HandleVisibility', 'off');

% Plot Initial Set reachable set 1 Plot Time = 1 Reachable Set RF(1, X0)
plot([max(x_T(1, :)), max(x_T(1, :))], [-1, 6], '--', ...
                            'Color', [0.6, 0.6, 0.6], ...
                            'LineWidth', 2, ...
                            'HandleVisibility', 'off');
plot([min(x_T(1, :)), min(x_T(1, :))], [-.1, +.1], 'Color', c, ...
                'LineWidth', 2, ... 
                'HandleVisibility', 'off');
plot([max(x_T(1, :)), max(x_T(1, :))], [-.1, +.1], 'Color', c, ...
                'LineWidth', 2, ... 
                'HandleVisibility', 'off');



for i = 1:2
    c = [0, 0, 1];
    SF = makeRectangle([xu_T1(:, i), xo_T1(:, i)]);
    plot([min(SF(1, :)), max(SF(1, :))], [i, i], 'Color', c, ...
                'LineWidth', 2, ... 
                'HandleVisibility', 'off');
    plot([min(SF(1, :)), min(SF(1, :))], [i-.1, i+.1], 'Color', c, ...
                'LineWidth', 2, ... 
                'HandleVisibility', 'off');
    plot([max(SF(1, :)), max(SF(1, :))], [i-.1, i+.1], 'Color', c, ...
                'LineWidth', 2, ... 
                'HandleVisibility', 'off');
end


Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

matlab2tikz('cos_reach1.tikz', 'width', '6cm', 'height', '4cm')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = d1(x, xh)
    out(1) = cos(x(2)) + x(2) - xh(2);

    out(2, 1) = 0; 
end

function out = d2(x, xh)
    out(1) = cos(xh(2)) + x(2) - xh(2);

    out(2, 1) = 0; 
end

function out = d3(x, xh)
    if all(x <= xh)
        out(1) = max([1, 0]*[d1(x, xh), d2(x, xh)]);
    elseif all(x >= xh)
        out(1) = min([1, 0]*[d1(x, xh), d2(x, xh)]);
    end

    out(2, 1) = 0; 
end

function out = d4(x, xh)
    out(1) = max(-1, min(1, [1, 0]*d3(x, xh)));
    out(2, 1) = 0; 
end

function out = d5(x, xh)
    out(1) = d_cos(x(2), xh(2));
    out(2, 1) = 0; 
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
    out = [ cos(x(2));...
            0 ];
end