

clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = 'cos';
name_latex = ''; %\cos


% dax2_bx_c(w, what, a);
f = @(w) d(w, w); 
tikz_option = 0;


%d = @(w, what) d_cosn(w, what, -3.1);
%f = @(w) cos(w)^(-3.1);

X_range = [-3, 3];




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



holder2 = [];
for i = 1:size(points, 2)
    xnow = points(:, i);
    out = d2(xnow(2), xnow(1));
    holder2 = [holder2, [xnow; out]];    
end
Xh2 = reshape(holder2(1, :), [num, num]);
X2 = reshape(holder2(2, :), [num, num]);
Z12 = reshape(holder2(3, :), [num, num]);



figure(1); clf;
hold on; grid on;
Leg = legend();
set(Leg,'visible','off');
%axis([X_range(1), X_range(2), ...
%      X_range(1), X_range(2), ...
%      (max(max(Z1)) - min(min(Z1)) +.001)*[-.05, .05]+[min(min(Z1)), max(max(Z1))]])
axis([-3, 3, -3, 3, -10, 10])
xticks([-3, 0, 3])
yticks([-3, 0, 3])
zticks([-10, -5, 0, 5, 10])
%if max(max(Z1)) == min(min(Z1))
%    zticks(min(min(Z1)))
%else
%    zticks(round([min(min(Z1)), (max(max(Z1)) + min(min(Z1)))/2, max(max(Z1))], 2))

%end
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
%view([.5, -1, 1])
%view([-.1, -1, .5])
view([48.0106, 28.2618])
xlabel('$\widehat{w}$','Interpreter','latex')
ylabel('$w$','Interpreter','latex')
zlabel(['$d^{\,', name_latex, '}(w, \widehat{w})$'],'Interpreter','latex')




% plot decomposition function 
s = surf(Xh2, X2, Z12, 'FaceColor','interp','edgecolor', 'none', 'HandleVisibility','off');

plot3(real(1, :), real(2, :),real(3, :), 'r', 'LineWidth', 2, 'HandleVisibility','off');

%plot3([-3, 3], [-3, 3], [-3, -3], 'k', 'LineWidth', 1.2, 'HandleVisibility','off');


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


drawnow

matlab2tikz('decomp2.tikz', 'width', '6cm', 'height', '4cm')


function out = d(x, xh)
    if x >= 0 && x >= - xh 
        out(1, 1) = x^2 + 2;   
    elseif xh <= 0 && x <= -xh
        out(1, 1) = xh^2 +  2;  
    elseif (x <= 0) && (xh >= 0)
        out(1, 1) = x*xh + 2;
    end
end


function out = d2(x, xh)
    if x >= 0 && x >= - xh 
        out(1, 1) = x^2 + 2;   
    elseif xh <= 0 && x <= -xh
        out(1, 1) = xh^2 +  2;  
    elseif (x <= 0) && (xh >= 0)
        out(1, 1) = 2;
    end
end



