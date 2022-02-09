%% -- Prior: rho(m) ---

function prior = prior_function(m,C,count)   %[rho,freq,n_events] = prior_function(m,C)   


 %Find patterns (configurations) for the upper left square of the model
n=0;
for k=1:16
    if  isequal(m(1:2,1:2),C(:,:,k))==1  
n=n+1;
k_initial(n)=k;
    end
end

freq=count(k_initial); %Compare that pattern with the training image distribution (hist)

%Conditional Probability of the chosen pattern
PcondA=count(k_initial)/sum(count);  %1/sum(freq); 


clear freq
%---FIRST 2 ROWS OF THE NEW IMAGE------
c=1; % If c=1, the 2x2 Matrix moves with overlap. Use always c=1!  
temp=0;
for i=1:c:1 %length(m(:,1))-1   
for j=1:c:length(m)-1
      
 temp=temp+1;
 
 %Find patterns (configurations) whose 1st column is the same as the current pattern
n=0;
for k=1:16
    if  isequal(m(i:i+1,j),C(:,1,k))==1  %isequal(C(:,1,RandV(temp)),C(:,1,k))==1 
n=n+1;
k_col(n)=k;
    end
end

freq=count(k_col); 

 %Find the 2x2 pattern which is exactly the same as the current square in the model
for n=1:length(k_col)
     if m(i:i+1,j:j+1)==C(:,:,k_col(n))
    k_equal=k_col(n);   
    end
end
 
 
    %Conditional Probability of the current pattern
    PcondB(temp)=count(k_equal)/sum(freq);
    
    
end
end



%%
%---REST OF THE NEW IMAGE-----
clear freq
clear k_equal
temp_step2=1;
temp_step3=1;

for i=2:c:length(m(:,1))-1   
for j=1:c:length(m)-1
        
  
   if j==1 %Starting a new row (left hand side)
       
  %--STEP2: Find patterns whose 1st ROW is the same as the current pattern---
n=0;
for k=1:16
    if  isequal(m(i,j:j+1),C(1,:,k))==1 
n=n+1;
k_row(n)=k;
    end
end

freq=count(k_row); 

  %Find the 2x2 pattern which is exactly the same as the current square in the model
for n=1:length(k_row)
     if m(i:i+1,j:j+1)==C(:,:,k_row(n))
    k_equal=k_row(n);   
    end
end

 
     %Conditional Probability of the chosen pattern
    PcondC(temp_step2)=count(k_equal)/sum(freq);

    temp_step2=temp_step2+1;
    
   else  %j>1  
   %--STEP3: Find 2x2 patterns whose unique unknown pixel is the botton rigth one---   
       n=0;
for k=1:16
    if  isequal(m(i,j:j+1),C(1,:,k))==1  &&  isequal(m(i+1,j),C(2,1,k))==1
n=n+1;
k_3(n)=k;
    end
end

freq3=count(k_3); 
     
 %Find the 2x2 pattern which is exactly the same as the current square in the model
for n=1:length(k_3)
     if m(i:i+1,j:j+1)==C(:,:,k_3(n))
    k_equal3=k_3(n);   
    end
end


      %Conditional Probability of the chosen pattern
    PcondD(temp_step3)=count(k_equal3)/sum(freq3);
    
     temp_step3=temp_step3+1;
    
   end 
   
  
    
end
end

%Conditional Probability of the entire model (prior)
prior=PcondA*prod(PcondB)*prod(PcondC)*prod(PcondD);


m0=m;   %m=m_new;



%{
%------Moving Window -> 2x2 Matrix ---------
freq=zeros(1,16);
n=0; %freq how many 2x2 patterns I compare in the original image

c=1; % If c=1, the 2x2 Matrix moves with overlap 

for i=1:c:length(A(:,1))-1   
for j=1:c:length(A)-1
    
    n=n+1;
    a=A(i:i+1,j:j+1);


 for k=1:16    %Compare the small matrix "a"2x2 with the matrix of all possible Configurations "C"
 
    if isequal(a,C(:,:,k))==1
 
 freq(1,k)=freq(1,k)+1; %freq when you find the pattern of a specific k-configuration
    end
        
 end


end
end

%Marginal Probability of a given pattern
 P=freq./n;

%--Pick Random configurations based on the Histogram (probability distribution)---
X=1:16; % possible values (configurations)
RandV=RandHistValues(freq,X); %Function

rho=freq(RandV(1))./n; % ??  I'm asuming than the prior of the whole new image (model) is a random P(x) that I select from the distribution (Histogram) of that model
n_events=n;
%}
end