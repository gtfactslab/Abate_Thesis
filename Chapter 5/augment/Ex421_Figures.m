

clc; clear all;

X0 = [0 , 1; ...
      0 , 1; ...
      0,  1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .0008;  % 
T  = .5;     % Simulation time 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1
figure(1); clf; hold on; grid on; 
axis([-.5 3 -.5 4 -1 4])
xlabel('$x_1$','Interpreter','latex')
%xticks([-6, -4, -2, 0, 2])
ylabel('$x_2$','Interpreter','latex')
%yticks([-2, 0, 2, 4])
zlabel('$x_3$','Interpreter','latex')
%zticks([-6, -4, -2, 0, 2])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
view([.5 -.5 .2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disc = 22;
holder1 = linspace(X0(1, 1), X0(1, 2), disc);
holder2 = linspace(X0(2, 1), X0(2, 2), disc);
holder3 = linspace(X0(2, 1), X0(2, 2), disc);
Start_Set = [];
for i =1:1:disc
    for j =1:1:disc
        for k =1:1:disc
            Start_Set = [Start_Set, [holder1(i); holder2(j); holder3(k)]];
        end
    end
end
clear holder1 holder2 holder3



holder = Start_Set;
holder2 = [];

Embedding_Traj = X0(:);
Embedding_Traj2 = [0; 0; 0; 0; 1; 1; 1; 2];


figure(1)
% Simulate system dynamics for reachable set
for t = 0:dt:T
    t
    if t ~= 0
        holder = Next_Set;
    end
    
    for i = 1:size(holder, 2)
        x_now = holder(:, i);
        x_next = x_now + dt*dxdt(x_now);
        holder2 = [holder2, x_next];
    end
    Next_Set = holder2;
    holder2 = [];

    Embedding_Traj = Embedding_Traj + dt*e(Embedding_Traj);
    Embedding_Traj2 = Embedding_Traj2 + dt*e2(Embedding_Traj2);
    
    if t ~= 0
        delete(a)
    end

    a(1) = scatter3(Next_Set(1, :), Next_Set(2, :), Next_Set(3, :), ...
                    'b', 'filled', 'HandleVisibility', 'off');

    drawnow
end
delete(a)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Next_Set_plot = Next_Set(:, 1:16:end);
k = boundary(Next_Set_plot', .8);
trisurf(k, Next_Set_plot(1, :)', Next_Set_plot(2, :)', Next_Set_plot(3, :)',...
        'FaceColor','green','FaceAlpha', 1, 'EdgeColor', [0, .6, 0], ...
        'HandleVisibility', 'off')

 Embedding_Traj = [Embedding_Traj(1:3), Embedding_Traj(4:6)];
 Embedding_Traj2 = [Embedding_Traj2(1:4), Embedding_Traj2(5:8)];

 holder = [];
 holder2 = [];
 for i = 1:2
     for j = 1:2
         for k = 1:2
             holder = [holder; [Embedding_Traj(1, i), Embedding_Traj(2, j), Embedding_Traj(3, k)]];
             holder2 = [holder2; [Embedding_Traj2(1, i), Embedding_Traj2(2, j), Embedding_Traj2(3, k)]];
          end
     end
 end
 k = boundary(holder, 0);
 trisurf(k, holder(:, 1), holder(:, 2), holder(:, 3), ...
         'FaceColor','blue','FaceAlpha',.1, 'HandleVisibility', 'off')
 trisurf(k, holder2(:, 1), holder2(:, 2), holder2(:, 3), ...
         'FaceColor','red','FaceAlpha',.1, 'HandleVisibility', 'off')

     
     
% Plot initial set
holder1 = linspace(X0(1, 1), X0(1, 2), 2);
holder2 = linspace(X0(2, 1), X0(2, 2), 2);
holder3 = linspace(X0(3, 1), X0(3, 2), 2);
Initial_Set = [];
for i =1:2
    for j =1:2
        for k =1:2
            Initial_Set = [Initial_Set; holder1(i), holder2(j), holder3(k)];
        end
    end
end
k_init = boundary(Initial_Set, 0);
trisurf(k_init ,Initial_Set(:, 1), Initial_Set(:, 2), Initial_Set(:, 3),...
    'FaceColor','red','FaceAlpha',1, 'HandleVisibility', 'off')


Leg = legend();
set(Leg,'visible','off')
drawnow

matlab2tikz('F5a.tikz', 'width', '6cm', 'height', '4cm')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compress to x1 x2 plane

figure(2); clf; hold on; grid on; 
axis([-.5 2 -.5 3])
xlabel('$x_1$','Interpreter','latex')
xticks([0, 1, 2])
ylabel('$x_3$','Interpreter','latex')
yticks([0, 1, 2, 3])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')


a = Embedding_Traj([1, 3], :);
over_approx = [a(1, 1), a(1, 2),a(1, 2),a(1, 1); ...
               a(2, 1), a(2, 1),a(2, 2),a(2, 2)];
patch(over_approx(1, :), over_approx(2, :), 'b', 'FaceAlpha', .2)

a = Embedding_Traj2([1, 3], :);
over_approx = [a(1, 1), a(1, 2),a(1, 2),a(1, 1); ...
               a(2, 1), a(2, 1),a(2, 2),a(2, 2)];
patch(over_approx(1, :), over_approx(2, :), 'r', 'FaceAlpha', .2)


Initial_Set_2D = [X0(1, 1), X0(1, 2), X0(1, 2), X0(1, 1); ...
                  X0(3, 1), X0(3, 1), X0(3, 2), X0(3, 2)];
patch(Initial_Set_2D(1, :), Initial_Set_2D(2, :), 'r')

Reachable_Set = Next_Set([1, 3], :);
k = boundary(Reachable_Set', .8);
Reachable_Set = Reachable_Set(:, k);
patch(Reachable_Set(1, :), Reachable_Set(2, :), 'g', 'FaceAlpha', .7)

Leg = legend();
set(Leg,'visible','off')
drawnow


%matlab2tikz('F5b.tikz', 'width', '6cm', 'height', '4cm')


function out = dxdt(x)
    out = [-x(1, 1)*x(2, 1) + x(3, 1) ; ...
           2*x(1, 1)*x(2, 1) ; ...
           x(1, 1) + x(2)];
end


function out = e(a)
    
    x = a(1:3, 1);
    xh = a(4:6, 1);
    
    out = [d(x, xh); ...
           d(xh, x)];
end

function out = e2(a)
    
    x = a(1:4, 1);
    xh = a(5:8, 1);
    
    out = [d2(x, xh); ...
           d2(xh, x)];
end

function out = d(x, xh)
    if x(1) <= 0
        out(1, 1) = -x(1)*x(2) + x(3);
    elseif x(1) >= 0
        out(1, 1) = -x(1)*xh(2) + x(3);
    end
    
    if x(2) >= 0
        out(2, 1) = 2*x(1)*x(2);
    elseif x(2) <= 0
        out(2, 1) = 2*xh(1)*x(2);
    end
    
    out(3, 1) = x(1) + x(2);
end


function out = d2(x, xh)
    if x(1) <= 0
        out(1, 1) = -x(1)*x(2) + x(3);
    elseif x(1) >= 0
        out(1, 1) = -x(1)*xh(2) + x(3);
    end
    
    if x(2) >= 0
        out(2, 1) = 2*x(1)*x(2);
    elseif x(2) <= 0
        out(2, 1) = 2*xh(1)*x(2);
    end
    
    out(3, 1) = x(4);
    
    if prod(x <= xh)
        out(4, 1) = x(3) + min([ x(1)*x(2), xh(1)*x(2),x(1)*xh(2), xh(1)*xh(2)]);
    elseif prod(xh <= x)
        out(4, 1) = x(3) + max([x(1)*x(2), xh(1)*x(2),x(1)*xh(2), xh(1)*xh(2)]);
    end
end
