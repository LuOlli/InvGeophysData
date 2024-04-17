%Code for solving Inversion of Geophysical data course exercise 3.
clear all
close all
g = 6.67e-11;%Nm^2/kg, gravitational constant
x = (0:40)'*1e3;
zG = [0.001 1:100]'*1e3;%depth for calculating G
z = [0.001 1:99]'*1e3;
G = nan*ones(length(x), length(z)-1);
for n = 1:length(zG)-1
    G(:,n)= 2*g*(atan2(zG(n+1),x)-atan2(zG(n),x));%Forward matrix (G)
end
d_true = 2*g*100*(atan2(15000,x)-atan2(5000,x)); %'true', calculated values of delta g
std = 2e-10;% I use a larger error value to better see the effects 
rng(147395104);% set a fixed seed so that the error stays the same between code executions
d = d_true+(std)*randn(length(d_true),1);%"measurements"
[U,S,V] = svd(G);%singular value decomposition
figure(1)
plot(diag(S),'k.')%the singular value spectrum
set(gca, 'YScale', 'log')
ylabel('Singular value')
figure(2)
% for i=1:10:length(V)
col = [1 4 8 12 16 20 26 32 40];
for i=1:9
    subplot(3,3,i)
    plot(z,V(:,col(i)));
    title(['Column ', num2str(col(i))])
end
xlabel('Depth (km)')
ylabel('V')
%% no noise
figure(3)
p1 = [1 3 5 7];%p values
for pp =1:4
Sp1= S(1:p1(pp),1:p1(pp));%truncated singular value matrix
Vp1 = V(:,1:p1(pp));%truncated matrix V
Up1 = U(:,1:p1(pp));%truncated matrix U
Gp1 = Up1*Sp1*Vp1';%recalculated G
m_p1 = Vp1*(Sp1\(Up1'*d_true));%inverted model parameters usinf truncated SVD
d_p1 = Gp1*m_p1;
subplot(2,2,pp)
hold on
plot(x,d_true,'r.')
plot(x,d_p1,'k--')
hold off
xlabel('Distance (m)')
ylabel('d or d_{true}')
title(['No noise in data, p = ',num2str(p1(pp))])
ylim([0 1e-8])
end
legend('Data','d_{p}')


%% added noise 
figure(4)
p2 = [3 5 10 20];
for pp =1:4
Sp2= S(1:p2(pp),1:p2(pp));%truncated singular value matrix
Vp2 = V(:,1:p2(pp));%truncated matrix V
Up2 = U(:,1:p2(pp));%truncated matrix U
Gp2 = Up2*Sp2*Vp2';%recalculated G
m_p2 = Vp2*(Sp2\(Up2'*d));%inverted model parameters usinf truncated SVD
d_p2 = Gp2*m_p2;
subplot(2,2,pp)
hold on
plot(x,d_true,'b')
plot(x,d,'r.')
plot(x,d_p2,'k--')
hold off
xlabel('Distance (m)')
ylabel('d or d_{true}')
title(['Noisy data, p = ',num2str(p2(pp))])
% ylim([0 1e-8])
end
legend('True data','Data','d_{p}')
%% model resolution matrices
p3 = [1 2 3 4 5 8 12 16 20 26 32 40];
figure(5)
hold on
for pp =1:length(p3)
Vp3 = V(:,1:p3(pp));
R = Vp3*Vp3';
plot(z,diag(R),'DisplayName',['p =',num2str(p3(pp))])
end
% set(gca, 'YScale', 'log')
legend('show')
hold off







