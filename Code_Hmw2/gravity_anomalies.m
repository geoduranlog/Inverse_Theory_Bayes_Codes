function g = gravity_anomalies(x_obs,x,m)

% INPUTS:
% x: horizontal discretization (vector)
% m : model parameters (thickness)

% OUTPUT:
% g : anomalie values at x


G_const = 6.674e-6; % Gravitational constant
delta_rho = -1700; % density contrast
deltax = x(2)-x(1); % grid size
delta = eps;

g = zeros(size(x_obs));

for i = 1: length(x_obs) 
    
    for j = 1:length(x)       
            
                  
            g(i) = g(i) + G_const * delta_rho * deltax * log(((x(j)-x_obs(i))^2 + m(j)^2)/((x(j)-x_obs(i))^2 + delta) );             
        
    end 
    
end




end
 