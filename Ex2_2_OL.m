%MATLAB code for solving exercise 2 problem 2.
g = 6.67e-11;%Nm^2/kg, gravitational constant
x = (0:40)'*1e3;
z = [0.001 1:100]'*1e3;
G = nan*ones(length(x), length(z)-1);
for n = 1:length(z)-1
    G(:,n)= atan2(z(n+1),x)-atan2(z(n),x);%Observation matrix (G)
end
d_true = 2*g*100*(atan2(15000,x)-atan2(5000,x)); %'true', calculated values of delta g
m_true = [zeros(1,5) 2*g*100*ones(1,10) zeros(1,85)]';% true values of m
std = 1e-8;%standard deviation for the added error
rng(5616518);% set a fixed seed so that the error stays the same between code executions
d = d_true+std*randn(length(d_true),1);%measurements with generated random Gaussian error
e=0.35;%dampening parameter
m_est =(G'*G+e*eye(length(z)-1))\(G'*d);%solving the estimated m with dampened LS
d_pre = G*m_est;%the predicted measurements
figure
subplot(1,2,1)
hold on
plot(m_true,z(2:end),'b')
plot(m_est,z(2:end),'r')
ylabel('Depth (m)')
xlabel('m or m_{true}')
xlim([-4e-9,16e-9])
ylim([0,1e5])
set(gca, 'YDir','reverse')
legend('m_{true}','m_{est}','Location','southeast')
hold off
subplot(1,2,2)
hold on
plot(x,d_true,'b')
plot(x,d_pre,'r--')
xlabel('Distance (m)')
ylabel('d or d_{true}')
ylim([-2e-9,12e-9])
xlim([0,4e4])
legend('d_{true}','d_{pre}','Location','southeast')
hold off