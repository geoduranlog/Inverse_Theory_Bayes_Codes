m_est=1.9e7;
m0=m_est;

m=(0.33:0.1:3)*m_est;

y=log(m/m0);

log(1/3);
sigma=log(3);

rho_y=exp( -( y.^ 2) / (2*sigma^2));

rho_m=m.*rho_y;

%Result 
figure(1)   
plot(y,rho_y)
set(gca,'fontsize',18)
grid on
title(['Probability Density Ln(m/mo) '],'fontsize',24);
xlabel('y')
ylabel('rho(y)')

%Result 
figure(2)   
plot(m,rho_m)
set(gca,'fontsize',18)
grid on
title(['Probability Density - mass '],'fontsize',24);
xlabel('m')
ylabel('rho(m)')
