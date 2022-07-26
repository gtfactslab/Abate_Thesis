
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 1/13/2021
% Description:  This script generates Figures 1a 1b and 1c.
%               Forward time reachable sets are predicted using MM.

clc; clear all;

global A Ap Am B Bp Bm
A = [1, 2; -1, 2];
Ap = [1, 2;  0, 2];
Am = [0, 0; -1, 0];
B = [0; 1];
Bp = [0; 1];
Bm = [0; 0];

W = [0, 1];


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
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r' , ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .1, ...
            'HandleVisibility', 'off');
% Plot MM Overapproximation of RF(1, X0) with d2
                  

% Plot Time = 1 Reachable Set RF(1, X0)
patch(x_T(1, :), x_T(2, :), 'g', ...
                            'FaceAlpha', .1, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');



Leg = legend();
set(Leg,'visible','off')

%matlab2tikz('emb1.tikz', 'width', '6cm', 'height', '4cm')


hi = [];
hr = [];
for j = 1:130
x10 = rand(2, 1) - .5;
X1 = [];
X1 = x10;
for i = 1:T_size
    xdot = dxdt(X1(:, i));
    X1(:, i+1) = X1(:, i) + dt*xdot;
end
hi = [hi, X1(:, 1)];
hr = [hr, X1(:, end)];
end
scatter(hi(1, :), hi(2, :), 15, [.8, 0, 0.1], 'filled', 'HandleVisibility', 'off')
scatter(hr(1, :), hr(2, :), 15, [.1, .6, 0.1], 'filled','HandleVisibility', 'off')
drawnow;

%matlab2tikz('reachable_sets_1.tikz', 'width', '4cm', 'height', '3cm')


figure(2); clf;
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
                  
ellipse(2,1,.6,2.8,1.3,'g',300)

% Plot Time = 1 Reachable Set RF(1, X0)
patch(x_T(1, :), x_T(2, :), [.1, .8, 0.1], ...
                            'FaceAlpha', .35, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');



Leg = legend();
set(Leg,'visible','off')



%matlab2tikz('reachable_sets_2.tikz', 'width', '3.5cm', 'height', '2cm')

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


function h=ellipse(ra,rb,ang,x0,y0,C,Nb)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semiminor axis of radius rb, and an orientation of the semimajor
% axis with an angle of ang (in radians) rotated counter-clockwise 
% from the x-axis.  The ellipse is centered at the point x0,y0.
%
% The length of ra, rb, and ang should be the same. 
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same 
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several ellipses C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be specified
% for each ellipse individually, in which case it should be the same
% length as ra, rb, etc.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% usage exmple: the following produces a red ellipse centered at 1,1
% and tipped down at a 45 deg axis from the x axis
% ellipse(1,2,pi/4,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%
% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original written by Peter Blattner, Institute of 
% Microtechnology, University of Neuchatel, Switzerland, blattner@imt.unine.ch
% Check the number of input arguments 
if nargin<1,
  ra=[];
end;
if nargin<2,
  rb=[];
end;
if nargin<3,
  ang=[];
end;
if nargin<5,
  x0=[];
  y0=[];
end;
 
if nargin<6,
  C=[];
end
if nargin<7,
  Nb=[];
end
% set up the default values
if isempty(ra),ra=1;end;
if isempty(rb),rb=1;end;
if isempty(ang),ang=0;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
if isempty(C),C=get(gca,'colororder');end;
% work on the variable sizes
x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);
if isstr(C),C=C(:);end;
if length(ra)~=length(rb),
  error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end;
% how many inscribed elllipses are plotted
if length(ra)~=length(x0)
  maxk=length(ra)*length(x0);
else
  maxk=length(ra);
end;
% drawing loop
for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    radm=ra(k);
    radn=rb(k);
    if length(ang)==1
      an=ang;
    else
      an=ang(k);
    end;
  elseif length(ra)==1
    xpos=x0(k);
    ypos=y0(k);
    radm=ra;
    radn=rb;
    an=ang;
  elseif length(x0)==length(ra)
    xpos=x0(k);
    ypos=y0(k);
    radm=ra(k);
    radn=rb(k);
    an=ang(k)
  else
    rada=ra(fix((k-1)/size(x0,1))+1);
    radb=rb(fix((k-1)/size(x0,1))+1);
    an=ang(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end;
  % draw ellipse
  
  co=cos(an);
  si=sin(an);
  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
  %  x=radm*cos(the)*co-si*radn*sin(the)+xpos;
  %  y=radm*cos(the)*si+co*radn*sin(the)+ypos;
  q1 = radm*cos(the)*co-si*radn*sin(the)+xpos;
  q2 = radm*cos(the)*si+co*radn*sin(the)+ypos;
  
  patch(q1, q2, C, 'FaceAlpha', .05, 'LineWidth', 1.25)
  
  % output handles to each ellipse if output variable specified
  
  if nargout > 0
    h(k)=p;
  end
  
end;
end

