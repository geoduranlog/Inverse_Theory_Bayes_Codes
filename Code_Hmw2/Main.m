%%% Assignment 2 -- Inverse theory 

clear all
clc
close all

%% INPUTS:

file2save = './TEST'; %'./Results_MCMC_10000samp_searchUnif65negative_stdprior300_burn0_17param.mat';

x_min = 0; % starting horizontal position
x_max = 3420; % end horizontal position
nx = 25; %17

% define our parameterization:
dx = (x_max-x_min)/(nx-1);
x = x_min:dx:x_min+dx*(nx-1);
x = x(:);

% Observed data
x_obs = [535; 749; 963; 1177; 1391; 1605; 1819; 2033; 2247; 2461; 2675; 2889];
Obs_Data = [-15; -24; -31.2; -36.8; -40.8; -42.7; -42.4; -40.9; -37.3; -31.5; -21.8; -12.8]; % Measured gravity anomalies at x positions (miligal)


% Mean model and standard deviation for the prior:

G_const = 6.674e-6; % Gravitational constant [miligal * m^2/kg]
delta_rho = -1700; % density contrast [kg/m^3]
% 
m_prior = (sum(Obs_Data/(2*pi*(delta_rho)*G_const))/length(Obs_Data))*ones(size(x)) ; % Bouguer plate [m]: ???preferred model???

m_prior(1) = 0; % anomalie zero in the boundaries
m_prior(end) = 0;

std_m = 250; % standard deviation parameters (prior)

%MCMC parameters:

n_samples=100000 %5000 %100000; %Total number of samples

burn_in = 900; %1000; %Lenght of the burn-in period
%% Starting model:

m=m_prior;

%Compute the likelihood and prior for the initial model

L = Likelihood(x,m,Obs_Data,x_obs);

prior = prior_function(m,m_prior,std_m);

%% MCMC code:

m_array=zeros(n_samples,length(m)+2);

n=1;
k_accepted=1;

for i_sample=1:n_samples    
      
    m_new = generate_model(m,50); %generate_model(m,15); % 50 %Set step size
       
    L_new = Likelihood(x,m_new,Obs_Data,x_obs);
    
    prior_new = prior_function(m_new,m_prior,std_m);
       
L_save(i_sample)=L; %Save the likelihood of each iteration %ADN

    
    if(rand<(L_new*prior_new)/(L*prior))
                
        m = m_new;
        L = L_new;
        prior = prior_new;
                
        posterior_accepted(k_accepted)=prior*L; %Save the posterior only for the accepted models %ADN
        
        k_accepted=k_accepted+1;
    end
    
    if (burn_in <= i_sample)
            m_array(n,:)=[m',L,prior]; %Each row  contains m_k (parameters of model m_k), its L and its prior. ADN
            n=n+1;
   end
    

end



[max,ind] = max(m_array(:,end-1).*m_array(:,end)); %Max of the posterior ADN

figure(2),plot(x,m_array(ind,1:length(x)),'r'),hold on,plot(x,m_prior,'k')
xlabel('x [m]')
ylabel('Thickness [m]')
axis ij

g = gravity_anomalies(x_obs,x,m_array(ind,1:length(x))); %g computed from the best found model (the one which maximize the posterior)
g_init = gravity_anomalies(x_obs,x,m_prior);

figure(3), plot(x_obs,g),hold on,plot(x_obs,Obs_Data) %,hold on, plot(x_obs,g_init,'--k')  %ADN
title('Data and Model-Data ')
set(gca,'fontsize',18)
grid on
xlabel('x (m)')
ylabel('dg [mGal]')
%legend('Synthetic','Observed Data' ,'Initial g')
legend('Synthetic Data','Observed Data')



%=======SAVE================
%save(file2save,'m_array')

%--Save Likelihoods vectors--
%L5=L_save;
%save L5 

%%
% L_plot=zeros(5,length(L_save));
% for i=1:5
% load (['L',num2str(i,'%01.0f'),'.mat']);   %Careful, this is loading all
% constants from when it was computed
%     
% %  L_plot(i,:)=(['L',num2str(i,'%01.0f')])  ;
% %  
% %  
% % %L_temp=load (['L',num2str(i,'%01.0f'),'.mat']);  
% % %L_plot(i,:)=L_temp;
% % 
% % L_plot(i,:)=load (['L',num2str(i,'%01.0f'),'.mat']);  
% % %clear L_temp
% 
% end

%L_plot=[L1;L2;L3;L4;L5];

% figure(4)
% for i=1:5
% plot(-log(L_plot(i,:)) )
% hold on
% end
% %plot(-log(L_save)) 
% set(gca,'fontsize',18)
% grid on
% title(['Convergence Properties '],'fontsize',24);
% xlabel('iterations')
% ylabel('-log(L)')
% %legend('step1','step2','step3','step4','step5')



%error = posterior -> for each model
e=m_array(ind,end-1).*m_array(ind,end); %Each h has a distribution, here I took the 
Avg_h=mean(m_array(:,1:end-2));  %Mean of each h after tried all models
%Avg_h=m_array(ind,1:end-2);   % Best model

%How many times the expected value is repeated? 
%How many times the values of my best model are repeated? 
figure(5)
plot(x,Avg_h,'-')
set(gca,'fontsize',18)
grid on
title(['Best model'],'fontsize',24);
xlabel('x (m)')
ylabel('h (m)')
legend('Best model','Data')


figure(6)
%hist(m_array(:,:))
hist(m_array(:,26),50)
set(gca,'fontsize',18)
grid on
title(['Histogram '],'fontsize',24);
xlabel('h (m)')
ylabel('frequency')
 

posterior=m_array(:,end-1).*m_array(:,end);
%posterior=posterior_accepted;

figure(7)
plot(posterior,'*')
%plot(m_array(:,end),m_array(:,end-1),'*')
%hist(posterior)
%hold on
%plot( mean(m_array(1:end-2,:)) )
%errorbar( Avg_h ,e)
%plot(Avg_h)
%hist(m_array(:,2))
set(gca,'fontsize',18)
grid on
title(['A posteriori distribution'],'fontsize',24);
xlabel('m')
ylabel('posterior')

N=n_samples;
figure(8);
scatterhist(m_array(end,:),m_array(1:end,11)); hold on; 
scatterhist(posterior,m_array(1:end,11)); hold on; 
%scatterhist(m_array(1:end,10),m_array(1:end,11),[ceil(sqrt(N)) ceil(sqrt(N))]); hold on; 
%contour(X,Y,Z,22,'--','LineWidth',2); colormap summer;
