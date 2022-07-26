
clc; clear all;

x = 0:.05:2;
for i = 1:size(x, 2)

    y(i) = 3*x(i)^2 - 6*x(i) + 4;
end

figure(1); clf;
hold on; grid on;
Leg = legend();
set(Leg,'visible','off');
axis([0, 2, 0, 4])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
xlabel('$p_{x}$','Interpreter','latex')
ylabel('$p_{y}$','Interpreter','latex')

plot(x, y, 'b', 'LineWidth', 2)

plot([.75, 1.25], [.1, .1], 'r', 'LineWidth', 2)

scatter(1,1, 80, 'filled', 'r')

matlab2tikz('parab1.tikz', 'width', '6cm', 'height', '4cm')


