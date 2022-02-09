function [m,Pcond] = generate_model2(C,count)    %(m_init,C,std)

%--Pick Random configurations based on the Histogram (probability distribution)---
X =1:16; 
RandV=RandHistValues(count,X); %Function

%N=sum(count); %N is the max # events (how many times I scanned the image with a 2x2 square)


%% --Making a New Image choosing random 2x2 patterns but based on the distribution---following homework order
m=zeros(41,41); %Size of my desired model

%---FIRST 2 ROWS OF THE NEW IMAGE------
%Set the upper left corner of the new image
m(1:2,1:2)=C(:,:,RandV(1));    %m(1:6,1:12)   

%Cond Prob. of 1st selection: Probability=count(RandV(1))/sum(count)=frequency of that pattern/frequency of all patterns
PcondA=count(RandV(1))/sum(count);

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

%----Pick Random configurations based on the subset of the Histogram -----
Randfreq=RandHistValues(freq,k_col);
 
 
    m(i:i+1,j+1)=C(:,2,Randfreq(temp)); 
    
    %Conditional Probability of the chosen pattern
    PcondB(temp)=count(Randfreq(temp))/sum(freq);
    
    
end
end


%%
%---REST OF THE NEW IMAGE-----
clear freq
clear Randfreq
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

%----Pick Random configurations based on the subset of the Histogram -----
Randfreq=RandHistValues(freq,k_col);
 

    m(i+1,j:j+1)=C(2,:,Randfreq(1));
 
     %Conditional Probability of the chosen pattern
    PcondC(temp_step2)=count(Randfreq(1))/sum(freq);

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
     
%----Pick Random configurations based on the subset of the Histogram -----
Randfreq3=RandHistValues(freq3,k_3);
 

    m(i+1,j+1)=C(2,2,Randfreq3(1));

      %Conditional Probability of the chosen pattern
    PcondD(temp_step3)=count(Randfreq3(1))/sum(freq3);
    
     temp_step3=temp_step3+1;
    
   end 
   
  
    
end
end

%Conditional Probability of the entire model
Pcond=PcondA*prod(PcondB)*prod(PcondC)*prod(PcondD);

end % End Function