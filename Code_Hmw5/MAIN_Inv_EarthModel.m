% Hmw 5 - Probabilistic Inverse THeory --
% J.A. Duran

clear all
close all

tic

%% --Load Data and training image--
data=load('seisdata.txt','-ascii');
pointspread=load('pointspread.txt','-ascii');

Im=imread('Strebelle.tif'); % TIFF extension is a multi-image file 
A=Im(:,:,3);  % The import image has size nxnx4 -> 1 to 3 show the image, and 4 seems to be only the background

%--Make a binary image--Optional (Use it for Australia image)
A = imbinarize(A);

% 
% figure(1)
% imagesc(data)   
% %imagesc(pointspread)   
% colorbar

%% ---Create matrix C 2x2x16 with all posible configurations---
%b=0;  w=255; %Black and White index number (according to the original image)
b=0;  w=1; %For binary image (Australia)

Per= permn([b w],2) ; % Matrix with posible rows configurations 

fil = permn([1:4],2) ; %Filas ->posible combinations of the rows made in the previous matrix

for k=1:16 %Final matrix. With all 16 options. Each one is a 2x2 matrix
C(:,:,k)=[Per(fil(k,1) ,:) ; Per(fil(k,2) ,:)];
end

%--Load Histogram info of the original image (see Hmw 4)------
%N is the max # events (how many times I scanned the image with a 2x2 square) N=sum(count)
%count is the frequency vector(how many times a given pattern C is repeated) 
load hist_Strebelle.mat  %load hist_Australia.mat  

figure(2)
bar(count/sum(count)) %plot(P,'*') %hist(P)
%title('Histogram - training Image','fontsize',14);
title('Marginal Probability distribution','fontsize',14); %('Probability of each possible configuration','fontsize',14);
xlabel('Configuration number')
ylabel('Marginal probability');  % ylabel('Frequency');
grid on
set(gca,'fontsize',16)
xlim([0 17])



%=========================================================
%% --- Initial model, prior & Likelihood  ---
 %Initial Model - Obtain a new model based on the training image
%m=generate_model(A,C,count) ;  
[m,prior_o] =generate_model2(C,count) ;
prior_o  %show initial prior it must be non zero (always)
mo=m; %jut to plot later the initial model

 %imshow(m)  
 %imagesc(m)
%=========================================================
%Standard deviation of noise is 7% of the largest data point
std=0.07*max(max(data));

syndata = conv2(double(m),pointspread,'same'); % conv2(A,pointspread,'same');

figure(2)   %figure(3) imagesc(data)  
imagesc(syndata)   
colorbar



%===LIKELIHOOD ===
%L=Likelihood(syndata,data,std) ; 
mis=misfit(syndata,data,std) ;  %To avoid rounds. Small values tends to 0 otherwise in L(m)
%=========================================================


 
 %% --- MONTECARLO Code  ---
 
 %---Parameters Montecarlo--
 n_samples=10000 %50000  %00; %How many iterations you run montecarlo
 burn_in = 2000; %400; %Check it in the Convergence image
 

n=1;
k_accepted=1;

for i_sample=1:n_samples    
     
   
     %Generate new model. 1st by selecting 1 pixel randomly and then,
     %choosing the opposite color (black or white) in that pixel           
   for l=1:100 %Loop to ensure prior_new is non zero
     m_new=m;
     IND=randi(numel(m_new));
     if m_new(IND)==1
     m_new(IND) =0; %round(rand); 
     else
      m_new(IND)=1;
     end
     
     %% Prior %-use it only if it is non zero
     
  %--Find coord of the random selected pixel--
    s=size(m); %Size of the matrix
    [fil,col] = ind2sub(s,IND); 
   

    %--PRIOR of the Perturbed Model----
    
    %CondProb of the neiborhood around the selected pixel -do it only if it is non zero
   %Condprob=Condprob_function(m,fil,col,C,count);
   %Condprob_new=Condprob_function(m_new,fil,col,C,count);  
   %prior_new=prior*Condprob_new/Condprob; %Avoid the zeros! 
   
   prior = prior_function(m,C,count); 
   prior_new = prior_function(m_new,C,count);

     if (prior_new~=0 && prior~=0)
        break
     end
   end
    %% Synthetic data
    syndata = conv2(double(m_new),pointspread,'same'); % conv2(A,pointspread,'same');
    
    %Likelihood
    %L_new=Likelihood(syndata,data,std);
    mis_new=misfit(syndata,data,std);
    
    
 
       
    %L_save(i_sample)=L_new; %Save the likelihood of each iteration 
    mis_save(i_sample)=mis_new; %Save the misfit of each iteration 

    
    %--Condition Metropolis---
   if(rand<= min(1,exp(mis_new-mis)*(prior_new/prior)) )
                
        m = m_new;
        mis = mis_new;%L = L_new;
        prior = prior_new;
                
       % post_accepted(k_accepted)=prior*exp(mis); %prior*L; %Save the posterior only for the accepted models 
        
        k_accepted=k_accepted+1;
    
    
    if (burn_in <= i_sample) % Take info only after burn_in region
   %Save each accepted model after the burn in iteration - Probably too big file but I need the models which maximize the posterior
        m_accepted(:,:,n)=m; 
       post_accepted(n)=prior*exp(mis);
       pseudopost_accepted(n)=prior*(-mis);
        % mis_accepted(n)=mis; %L_accepted(n)=L;
       % prior_accepted(n)=prior;
            n=n+1;
   end
  end

end

%Max of the posterior - Then, the best model is m_accepted(:,:,ind)
 %[max,ind] = max(L_accepted.*prior_accepted); 
[max,ind] = max(post_accepted);
[min,ind2] = min(pseudopost_accepted);

%==SAVE THE MODELS==
%save('m_accepted_it50mil','m_accepted','mo','mis_save','ind')

%Convergence analysis
figure(4)   
plot(-mis_save);
set(gca,'fontsize',18)
grid on
title(['Convergence Properties '],'fontsize',24);
xlabel('iterations')
ylabel('|d-dpred|')




%--Result---
m_avg=mean(m_accepted,3);
figure(5)    
%imagesc(conv2(m_avg,pointspread,'same'));  
%imagesc(conv2(m_accepted(:,:,ind),pointspread,'same'));  
%imagesc(m_accepted(:,:,ind));  
imagesc(m_avg); 
colorbar


%--True Model - Objective---
Obj=load('truemodel.txt','-ascii');
%data2=conv2(Obj,pointspread,'same');
figure(6)
imagesc(data);
%imagesc(Obj)
colorbar


toc
 