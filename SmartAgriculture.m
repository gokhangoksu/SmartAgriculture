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
subplot(1,2,1);
plot(t,x);
title('Açýk Çevrim Sistemi')
xlabel('Zaman (saat)')
ylabel('Nem farký (mm)')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

subplot(1,2,2);
plot(t,xCLS);
hold on
title('Kapalý Çevrim Sistemi')
xlabel('Zaman (saat)')
ylabel('Nem farký (mm)')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

figure
u=-k*x;
hold on
plot(t,u);
plot(t,30*ones(size(t)),'r--');
plot(t,-30*ones(size(t)),'r--');
hold off
title('Kontrol Fonksiyonu')
xlabel('Zaman (saat)')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

figure
plot(t,w);
title('Bozucu Girdileri')
legend('R(t)','E(t)')
xlabel('Zaman (saat)')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
