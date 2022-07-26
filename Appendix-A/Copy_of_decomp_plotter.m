

clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [1, 0, 0];

name = 'cos';
name_latex = 'a \theta^2 + b \theta + c';
%d = @(w, what) (w<=what)*min([cos(w) + (w - what), cos(what) + (w - what)]) + ...
%(w> what)*min([cos(w) + (w - what), cos(what) + (w - what)]); 




d = @(w, what) dax2_bx_c(w, what, a);
% dax2_bx_c(w, what, a);
f = @(w) a(1)*w^2 + a(2)*w + a(3);
tikz_option = 0;


%d = @(w, what) d_cosn(w, what, -3.1);
%f = @(w) cos(w)^(-3.1);

X_range = [-2, 2];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE AND PLOT DECOMPOSTIION FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx1 = (X_range(2) - X_range(1))/(200); % diagnol plot discritization
dx2 = (X_range(2) - X_range(1))/(200); % surface plot discritization

[Xh, X] = meshgrid(X_range(1):dx2:X_range(2), ...
                   X_range(1):dx2:X_range(2));
num = size(Xh, 1);
               
points = [Xh(:)'; X(:)'];


real = [];
% get diagnol
for x = X_range(1):dx1:X_range(2)
    out = f(x);
    real = [real, [x; x; out]];
end

holder = [];
for i = 1:size(points, 2)
    xnow = points(:, i);
    out = d(xnow(2), xnow(1));
    holder = [holder, [xnow; out]];  
    %holder(3, i) = max([min([holder(3, i), 5.1])], -5.1);
end
Xh = reshape(holder(1, :), [num, num]);
X = reshape(holder(2, :), [num, num]);
Z1 = reshape(holder(3, :), [num, num]);

figure(1); clf;
hold on; grid on;
Leg = legend();
set(Leg,'visible','off');
%axis([X_range(1), X_range(2), ...
%      X_range(1), X_range(2), ...
%      (max(max(Z1)) - min(min(Z1)) +.001)*[-.05, .05]+[min(min(Z1)), max(max(Z1))]])
axis([-2, 2, -2, 2, 0, 4])
zticks([-5, 0, 5])
%if max(max(Z1)) == min(min(Z1))
%    zticks(min(min(Z1)))
%else
%    zticks(round([min(min(Z1)), (max(max(Z1)) + min(min(Z1)))/2, max(max(Z1))], 2))

%end
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
view([.5, -1, 1])
%view([-.1, -1, .5])
%view([1, .1, .5])
xlabel('$\widehat{\theta}$','Interpreter','latex')
ylabel('$\theta$','Interpreter','latex')
zlabel(['$d^{\,', name_latex, '}(\theta, \widehat{\theta})$'],'Interpreter','latex')

% plot decomposition function 
s = surf(Xh, X, Z1, 'FaceColor','interp', 'edgecolor', 'none', 'HandleVisibility','off');

plot3(real(1, :), real(2, :),real(3, :), 'r', 'LineWidth', 2, 'HandleVisibility','off');

%

%%Extract X,Y and Z data from surface plot
x=s.XData;
y=s.YData;
z=s.ZData;
%%Create vectors out of surface's XData and YData
x=x(1,:);
y=y(:,1);
%%Divide the lengths by the number of lines needed
xnumlines = 18; % 10 lines
ynumlines = 18; % 10 partitions
xspacing = round((length(x))/xnumlines);
yspacing = round((length(y))/ynumlines);
%%Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x, Y1,Z1,'-k', 'HandleVisibility','off');
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2, y,Z2,'-k', 'HandleVisibility','off');
end


plot3(real(2, :),real(1, :),real(3, :), 'r', 'LineWidth', 2, 'HandleVisibility','off');
drawnow

if tikz_option == 1
    matlab2tikz(['d_', name,'.tikz'], 'width', '6cm', 'height', '4cm')
end

