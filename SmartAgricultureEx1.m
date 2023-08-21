clear all
close all
clc

a=0.1;
c1=1.5;
c2=0.7;
C=[c1 -c2];

setlmis([]) 
p = lmivar(1,[1 1]);    % 1x1 sym
y = lmivar(1,[1 1]);    % 1x1 
gamma=0.94; % Feasible çýkan en küçük deðere kadar bulmaya çalýþ!!!

lmiterm([-1 1 1 p],1,1);        % 0 < p
lmiterm([2 1 1 p],-1,a,'s');    % -2ap < 0
lmiterm([2 1 1 y],-1,1,'s');    % -2ap-2y < 0
lmiterm([2 1 2 0],C);           % C < 0
lmiterm([2 1 3 p],1,1);         % p < 0
lmiterm([2 2 2 0],-gamma^2);    % -gamma^2*I < 0
lmiterm([2 3 3 0],-1);          % -1 < 0
lmis = getlmis;

[tmin,xfeas] = feasp(lmis);

p = dec2mat(lmis,xfeas,p);
y = dec2mat(lmis,xfeas,y);

k=y/p;

t=0:0.01:10;
w=[exp(-(t-2).^2)+exp(-(t-6).^2);
   ((t-4).^2+1).^(-1)+((t-8).^2+1).^(-1)];
x0=[10;8;6;4;2;-2;-4;-6;-8;-10];
x=zeros(length(x0),length(t));
xCLS=x;
for i=1:length(t)
    x(:,i)=exp(-a*t(i))*x0+C*w(:,i);
    xCLS(:,i)=exp(-(a+k)*t(i))*x0+C*w(:,i);
end

figure
subplot(1,19,1:8);
plot(t,x,'LineWidth',3.0);
title('Açýk Çevrim Sistemi')
xlabel('Zaman (saat)')
ylabel('Nem farký (mm)')
set(gca,'fontsize',18)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

subplot(1,19,11:18);
plot(t,xCLS,'LineWidth',3.0);
hold on
title('Kapalý Çevrim Sistemi')
xlabel('Zaman (saat)')
ylabel('Nem farký (mm)')
set(gca,'fontsize',18)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

ax = subplot(1,19,19,'Visible','off');
axPos = ax.Position;
delete(ax)
hL = legend('S_1','S_2','S_3','S_4','S_5','S_6','S_7','S_8','S_9','S_{10}');
hL.Position(1:2) = axPos(1:2)+[-0.02 0];
set(gca,'fontsize',18)

figure
u=-k*x;
subplot(1,10,1:9)
hold on
plot(t,u,'LineWidth',3.0);
plot(t,30*ones(size(t)),'r--','LineWidth',3.0);
plot(t,-30*ones(size(t)),'r--','LineWidth',3.0);
hold off
title('Kontrol Fonksiyonu')
xlabel('Zaman (saat)')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

ax = subplot(1,10,10,'Visible','off');
axPos = ax.Position;
delete(ax)
hL = legend('S_1','S_2','S_3','S_4','S_5','S_6','S_7','S_8','S_9','S_{10}');
hL.Position(1:2) = axPos(1:2)+[-0.05 -0.4];
set(gca,'fontsize',18)

figure
plot(t,w,'LineWidth',3.0);
title('Bozucu Girdileri')
legend('R(t)','E(t)')
xlabel('Zaman (saat)')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);