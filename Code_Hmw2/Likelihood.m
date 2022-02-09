%% Likelihood

function L = Likelihood(x,m,Obs_Data,x_obs)


synthetic = gravity_anomalies(x_obs,x,m);
sigma = 1; %standard deviation of the noise is 1.0 mgal

L = exp(-(norm( Obs_Data - synthetic ) ^ 2) / (2*sigma^2));


end