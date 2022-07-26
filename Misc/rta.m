




figure(1); clf; hold on;
axis([0, 10, 3, 7])


patch(-.5 + [1, 3, 3, 1], [5, 5, 6, 6], 'b', 'FaceAlpha', .05)

patch(-.5 + 3.5 + [1, 3, 3, 1], [4.5, 4.5, 6.5, 6.5], 'b', 'FaceAlpha', .05)

patch(-.5 + 7 + [1, 3, 3, 1], [5, 5, 6, 6], 'b', 'FaceAlpha', .05)

patch(-.5 + 3.5 + [1.15, 2.85, 2.85, 1.15], -.1+[4.75, 4.75, 5.5, 5.5], 'b', 'FaceAlpha', .05)


p1 = [2.5 5.5];                         % First Point
p2 = [4 5.5];                         % Second Point
dp = p2-p1;                         % Difference

quiver(p1(1),p1(2),dp(1),dp(2),1, 'Color','k','LineWidth',1.5)

p1 = [3.5+2.5, 5.5];                         % First Point
p2 = [3.5+4 5.5];                         % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1, 'Color','k','LineWidth',1.5)

p1 = [8.5, 5];                         % First Point
p2 = [8.5, 3.5];                         % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1, 'Color','k','LineWidth',1.5)

p1 = [1.5, 3.5];                         % First Point
p2 = [1.5, 5];                         % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1, 'Color','k','LineWidth',1.5)

p1 = [8.5, 3.5];                         % Second Point
p2 = [1.5, 3.5];                         % First Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1, 'Color','k','LineWidth',1.5)


text(.55,5.65,'Potentially Unsafe', 'FontSize',14)
text(.95,5.35,'Controller', 'FontSize',14)
text(4.45,5.15,'Prediction', 'FontSize',14)
text(4.4,4.9,'Mechanism', 'FontSize',14)
text(4.75,5.9,'RTA', 'FontSize',14)

text(8.2, 5.5,'Plant', 'FontSize',14)

matlab2tikz('rta_thing.tikz', 'width', '6cm', 'height', '4cm')

