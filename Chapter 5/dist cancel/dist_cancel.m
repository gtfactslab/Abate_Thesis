
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 9/9/2021
% Description:  This script generates Figures 1a and 1b.
%               Forward time reachable sets are predicted using MM.

clc; clear all; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Predict Reachable Sets using d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intervals Defining Initial Set
x0 = [1; 1];
W = [0, 1];
W_set1 = W(1):.5:W(2);
W_set2 = [0, .5, 1, 0, .5, 1, 0, .5, 1;...
          0, 0, 0, .5, .5, .5, 1, 1, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .005;   % Timestep for Simulation
T  = 1;      % Prediction Time-Horizon
T_size = size(0:dt:T, 2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on;
ax = gca;
%axis([-1, 5, -1, 3])
%xticks([-1, 0, 1, 2, 3, 4, 5])
%yticks([-1, 0, 1, 2, 3])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Initial Set X0
scatter(x0(1), x0(2), 80, 'r', 'filled')
drawnow

%%%%%%%%%%%%%%%%%%%%%%%
% Plot MM overapproximation
%%%%%%%%%%%%%%%%%%%%%%%
Tr =  [1, -1; 0, 1]; 

a = [x0; x0];
aT = [inv(Tr)*x0; inv(Tr)*x0];
for t = 1:T_size
    a = a + dt* [d(a(1:2, 1), W(1), a(3:4, 1), W(2)) ; ...
                 d(a(3:4, 1), W(2), a(1:2, 1), W(1))];
    aT = aT + dt*[dT(aT(1:2, 1), W(1), aT(3:4, 1), W(2)); ...
                  dT(aT(3:4, 1), W(2), aT(1:2, 1), W(1))];
end

% Plot MM Overapproximation of RF(1, X0)
SF = makeRectangle([a(1:2, 1), a(3:4, 1)]);
patch(SF(1, :), SF(2, :), 'g', ...
                          'FaceAlpha', .1, ...
                          'LineWidth', 1.25, ...
                          'HandleVisibility', 'off');
scatter([a(1), a(3)], [a(2), a(4)], 'k', 'filled', 'HandleVisibility', 'off');

SF = Tr*makeRectangle([aT(1:2, 1), aT(3:4, 1)]);
patch(SF(1, :), SF(2, :), 'r', ...
                          'FaceAlpha', .1, ...
                          'LineWidth', 1.25, ...
                          'HandleVisibility', 'off');

drawnow







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Reachable Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


holder = x0;
 
% Compute Time = 1 Second Reachable Set of System
% Compute MM approximation of Reachable Set
for t = 1:T_size
    t
    % reachable set computation
    holder2 = [];
    for i = 1:size(holder, 2)
        for j = 1:size(W_set2, 2)
            w_now = W_set2(:, j);
            x = holder(:, i);
            x_next = x + dt*dxdt2(x, w_now);
            holder2(:, end + 1) = x_next;
        end
    end

    if t >= 2
        k = boundary(holder2(1, :)', holder2(2, :)');
        holder = holder2(:, k);
    else
        holder = holder2;
    end 
end
% parse
holder5 = holder(:, 1);
for i = 2:size(holder, 2)
    if norm(holder5(:, end) - holder(:, i)) >= .05
        holder5(:, end + 1) = holder(:, i);
    end
end
x_T = holder5;

% Plot Time = 1 Reachable Set RF(1, X0)
patch(x_T(1, :), x_T(2, :), 'b', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');
drawnow





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Reachable Set Real
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

holder = x0;
 
% Compute Time = 1 Second Reachable Set of System
% Compute MM approximation of Reachable Set
for t = 1:T_size
    t
    % reachable set computation
    holder2 = [];
    for i = 1:size(holder, 2)
        for j = 1:size(W_set1, 2)
            w_now = W_set1(j);
            x = holder(:, i);
            x_next = x + dt*dxdt(x, w_now);
            holder2(:, end + 1) = x_next;
        end
    end

    if t >= 2
        k = boundary(holder2(1, :)', holder2(2, :)');
        holder = holder2(:, k);
    else
        holder = holder2;
    end 
end
holder5 = holder(:, 1);
for i = 2:size(holder, 2)
    if norm(holder5(:, end) - holder(:, i)) >= .02
        holder5(:, end + 1) = holder(:, i);
    end
end
x_T2 = holder5;

% Plot Time = 1 Reachable Set RF(1, X0)
patch(x_T2(1, :), x_T2(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');
drawnow

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%%
matlab2tikz('dist_cancel.tikz', 'width', '6cm', 'height', '4cm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = d(x, w, xh, wh)
    if x(2, 1) >= 0 && x(2, 1) >= - xh(2, 1)  
        out(1, 1) = x(2, 1)^2 + w;   
    elseif xh(2, 1) <= 0 && x(2, 1) <= -xh(2, 1)
        out(1, 1) = xh(2, 1)^2 +  w;  
    elseif (x(2, 1) <= 0) && (xh(2, 1) >= 0)
        out(1, 1) = x(2, 1)*xh(2, 1) + w;
    end
    
    out(2, 1) = x(1) - wh; 
end


function out = makeRectangle(X0)
    out = [X0(1, 1), X0(1, 2), X0(1, 2), X0(1, 1); ...
           X0(2, 1), X0(2, 1), X0(2, 2), X0(2, 2)];
end

function out = dxdt(x, w)
    out = [ x(2, 1)^2 + w;...
            x(1, 1) - w];
end

function out = dxdt2(x, w)
    out = [ x(2, 1)^2 + w(1);...
            x(1, 1) - w(2)];
end

function out = dT(x, w, xh, wh)
    % x2^2 - x2 + x1
    % x1 - w - x2
    
    out(1, 1) = dax2_bx_c(x(2), xh(2), [1; -1; x(1)]);
    out(2, 1) = x(1) - x(2) - wh; 
end




