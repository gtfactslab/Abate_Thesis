%from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000739



global delta k1 k2 eT

delta=1;
k1=1;
k2=1;
eT=1;

T=[0,100];

x0= [1.5,1.5,0,0];
x0(3:4)=[cos(0),sin(0)];

[t,xs]=ode45(@sys,T,x0);


%plot(t,xs(:,2));
figure(1); clf
grid on; hold on;
plot(xs(:,1),xs(:,2), 'b'); 








emb0=[x0+[-.5,-.5,0,0],x0+[0.5,0.5,0,0]];

emb0(3:4)=[cos(0),sin(0)];
emb0(7:8)=[cos(0),sin(0)];
wlow=[-.2,-.2];
whigh=[0.1,0.1];

% emb_aug0=[x0-[.01,.01],inf,x0+[0.01,0.01],inf];
% emb_aug0(3)=emb_aug0(1)*emb_aug0(2);
% emb_aug0(6)=emb_aug0(4)*emb_aug0(5);
% 
[te,embs]=ode45(@(t,y) [decomp(y(1:4),y(5:8),wlow,whigh); decomp(y(5:8),y(1:4),whigh,wlow)],T,emb0);
% [tea,emb_augs]=ode45(@(t,y) [decomp_aug(y(1:3),y(4:6)); decomp_aug(y(4:6),y(1:3))],T,emb_aug0);
% 
% 
figure(2); clf;
ax=subplot(2,1,1);
hold on; grid on;
fill([te;flipud(te)],[embs(:,1);flipud(embs(:,5))],'g','FaceAlpha',.2);
plot(t,xs(:,1));
ylabel('x');
xlabel('t')
ax=subplot(2,1,2);
hold on; grid on;
fill([te;flipud(te)],[embs(:,2);flipud(embs(:,6))],'r','FaceAlpha',.2);
plot(t,xs(:,2));
xlabel('t')
ylabel('y');
% 
% fill([tea;flipud(tea)],[emb_augs(:,1);flipud(emb_augs(:,4))],'g','FaceAlpha',.6);
% fill([tea;flipud(tea)],[emb_augs(:,2);flipud(emb_augs(:,5))],'r','FaceAlpha',.6);
% 
% 



function xdot=sys(t,x)
global delta k1 k2 eT
    xdot=[0;0];
    xdot(1)=x(3)-delta*x(1)+k1*x(2)-k2*(eT-x(2))*x(1);
    xdot(2)=-k1*x(2)+k2*(eT-x(2))*x(1);
    xdot(3)=-x(4);
    xdot(4)=x(3);
end

function out = decomp(x,xhat,w,what)
    global delta k1 k2 eT
    out=[0;0;0;0];
    out(1)=x(3)-delta*x(1)+k1*x(2)-k2*(eT-x(2))*x(1)+w(1);
    out(2)=-k1*x(2)+k2*(eT-x(2))*x(1)+w(2);
    out(3)=-xhat(4);
    out(4)=x(3);
end


    