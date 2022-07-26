
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
hold on; 
ax = gca;
axis([-1, 5, -1, 3])
xticks([-1, 0, 1, 2, 3, 4, 5])
yticks([-1, 0, 1, 2, 3])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r', ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .15, ...
            'HandleVisibility', 'off');
% Plot MM Overapproximation of RF(1, X0) with d2


% Plot Time = 1 Reachable Set RF(1, X0)
patch(x_T(1, :), x_T(2, :), [.1, .8, 0.1], ...
                            'FaceAlpha', .35, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');



Leg = legend();
set(Leg,'visible','off')

xo = [.5; .5];
Xo = zeros(2, T_size);
Xo(:, 1) = xo;
xu = -[.5; .5];
Xu = zeros(2, T_size);
Xu(:, 1) = xu;
tic
for i = 1:T_size
    Xu(:, i+1) = Xu(:, i) +dt*dxdt(Xu(:, i));
    Xo(:, i+1) = Xo(:, i) +dt*dxdt(Xo(:, i));
end
toc

xu_T = Xu(:, end);
xo_T = Xo(:, end);

rect = [xu_T(1, 1), xo_T(1, 1), xo_T(1, 1), xu_T(1, 1); ...
        xu_T(2, 1), xu_T(2, 1), xo_T(2, 1), xo_T(2, 1)]
% Plot Initial Set X0
patch(rect(1, :), rect(2, :), 'r' , ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .1, ...
            'HandleVisibility', 'off');

plot(Xu(1, :), Xu(2, :), 'b')
plot(Xo(1, :), Xo(2, :), 'b')
%matlab2tikz('nonmonotone.tikz', 'width', '3.5cm', 'height', '2cm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



