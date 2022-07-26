

clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [1, 1, -1];

name = 'cos';
name_latex = ''; %\cos
%d = @(w, what) F(w, what, a); 
%
d = @(w, what) d_cos(w, what); 




% dax2_bx_c(w, what, a);
f = @(w) d(w, w); %F(w, w, a); %a(1)*w^2 + a(2)*w + a(3);
tikz_option = 1;


%d = @(w, what) d_cosn(w, what, -3.1);
%f = @(w) cos(w)^(-3.1);

X_range = [-10, 10];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE AND PLOT DECOMPOSTIION FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx1 = (X_range(2) - X_range(1))/(50); % diagnol plot discritization
dx2 = (X_range(2) - X_range(1))/(20); % surface plot discritization

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
%axis([-10, 10, -10, 10, -5, 5])
zticks([-5, 0, 5])
%if max(max(Z1)) == min(min(Z1))
%    zticks(min(min(Z1)))
%else
%    zticks(round([min(min(Z1)), (max(max(Z1)) + min(min(Z1)))/2, max(max(Z1))], 2))

%end
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
%view([.5, -1, 1])
view([-.1, -1, .5])
view([1, .1, .5])
xlabel('$\widehat{w}$','Interpreter','latex')
ylabel('$w$','Interpreter','latex')
zlabel(['$d^{\,', name_latex, '}(w, \widehat{w})$'],'Interpreter','latex')

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
xnumlines = 10; % 10 lines
ynumlines = 5; % 10 partitions
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

view([48.0106,  28.2618])

drawnow

if tikz_option == 1
    matlab2tikz(['d_', name,'.tikz'], 'width', '6cm', 'height', '4cm')
end

function out = F(x, xh, v)
    a = v(1);
    b = v(2);
    c = v(3);
    if a*x <= -b/(2) && -b/(2) <= a*xh
        out = (-b^2 +4*a*c)/(4*a);
    end
    if a*x + a*xh + b <= 0 && ~(a*x <= -b/(2) && -b/(2) <= a*xh)
        out = a*xh^2 + b*xh + c; 
    end
    if 0 <= b + a*xh + a*x && ~(a*x <= -b/(2) && -b/(2) <= a*xh)
        out = a*x^2 + b*x + c;
    end
end
