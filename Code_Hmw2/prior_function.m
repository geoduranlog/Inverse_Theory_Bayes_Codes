%% Compute prior:

function prior = prior_function(m,m_prior,std)

  prior = exp(-0.5*norm(m-m_prior)^2/(std*std)); % new model, from uniform distribution with mean m_init and standard deviation std

end