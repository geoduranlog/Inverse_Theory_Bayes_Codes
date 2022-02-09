%%% Generate model from uniform distribution centered in the previous model

function m = generate_model(m_init,std)

std_vec = std.*ones(size(m_init));

std_vec(1) = 0; % fix the boundaries
std_vec(end) = 0;

% ind = find(m_init < std);
% std_vec(ind) = m_init(ind);

% m = m_init + std_vec.*rand(size(m_init)); % new model, from normal distribution with mean m_init and standard deviation std_vec

a = m_init - std_vec;
b = m_init + std_vec;


m = a + (b-a).*rand(size(m_init)); % new model, from uniform distribution with mean m_init and serch radious defined by std_vec

while any(m<0)
    m = a + (b-a).*rand(size(m_init)); % new model, from uniform distribution with mean m_init and serch radious defined by std_vec
end

end