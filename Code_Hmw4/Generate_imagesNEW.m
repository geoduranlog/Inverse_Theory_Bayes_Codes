% Hmw 4 - Probabilistic Inverse THeory --
% J.A. Duran

%{ 
Objective: Make an algorythn that randomly generates new images that are (statistically)
similar to a given training/reference image (of the Earth). Then this algorithm can later be used to
sample prior models in a Monte Carlo solution to an inverse problem.
%}

clear all
close all

%---Read and plot (save) the image------
%imfinfo('Strebelle.tif')
Im=imread('Australia.tif');
%Im=imread('Strebelle.tif'); % TIFF extension is a multi-image file 
%A=imread('Test3_Str.png'); %Only for computing the Marg Probability of the new images

A=Im(:,:,3);  % The import image has size nxnx4 -> 1 to 3 show the image, and 4 seems to be only the background

%figure(1)
%image(A)    % Here white look yellow and black look blue

figure(2)
imshow(A)    %Show the black- white image

%--Test if you can write it in an "easier" format
%imwrite(A,'Strebelle_simple.png') 


%--Make a binary image--Optional (Use it for Australia image)
BW = imbinarize(A);

%figure
%imshowpair(A,BW,'montage')

A=BW;   %A(170:178,208:219)
%image(A)

%%
%--Create matrix C 2x2x16 with all posible configurations---
%b=0;  w=255; %Black and White index number (according to the original image)
b=0;  w=1; %For binary image (Australia)

Per= permn([b w],2) ; % Matrix with posible rows configurations 

fil = permn([1:4],2) ; %Filas ->posible combinations of the rows made in the previous matrix

for k=1:16 %Final matrix. With all 16 options. Each one is a 2x2 matrix
C(:,:,k)=[Per(fil(k,1) ,:) ; Per(fil(k,2) ,:)];
end



%------Moving Window -> 2x2 Matrix ---------
count=zeros(1,16);
N=0; %Count how many 2x2 patterns I compare in the original image

c=1; % If c=2, the 2x2 Matrix moves without overlap. If c=1 there is overlap  

for i=1:c:length(A(:,1))-1   
for j=1:c:length(A)-1
    
    N=N+1;
    a=A(i:i+1,j:j+1);


 for k=1:16    %Compare the small matrix "a"2x2 with the matrix of all possible Configurations "C"
 
    if isequal(a,C(:,:,k))==1
 
 count(1,k)=count(1,k)+1; %Count when you find the pattern of an specific k-configuration
    end
        
 end


end
end

%======Probability of each configuration======
%N_nointer=( (length(A))^2 -1)/4; %Max Number of possible events for a 2x2 given configuration, with no overlapp betwwn patterns 
 P=count./N;  %OJO, if two values of count are too close => count./N could give the same result (Matlab rounds the numbers), better choose a most likely value based on the array "count" instead of P. Example P(14) and P(15) 

 %==SAVE HIST INFO - Training Image====
%save('hist_Strebelle.mat','count','N');   save('hist_Australia.mat','count','N'); 
 
figure(3)
bar(P) %plot(P,'*') %hist(P)
%title('Histogram - training Image','fontsize',14);
title('Marginal Probability distribution','fontsize',14); %('Probability of each possible configuration','fontsize',14);
xlabel('Configuration number')
ylabel('Marginal probability');  % ylabel('Frequency');
grid on
set(gca,'fontsize',16)
xlim([0 17])


%--Pick Random configurations based on the Histogram (probability distribution)---
% data obtained by hist -Configuration Number
X =1:16; % [.1 .4 .6 .8] ; % possible values (configurations)
%count % how many times they occurred

RandV=RandHistValues(count,X); %Function


%figure(4)
 %hist(RandV)

%% --Making a New Image choosing random 2x2 patterns but based on the distribution---following homework order


B=zeros(size(A));

%---FIRST 2 ROWS OF THE NEW IMAGE------
%Set the upper left corner of the new image
B(1:2,1:2)=C(:,:,RandV(1));    %B(1:6,1:12)

temp=0;
for i=1:c:1 %length(B(:,1))-1   
for j=1:c:length(B)-1
      
 temp=temp+1;
 
 %Find patterns (configurations) whose 1st column is the same as the current pattern
n=0;
for k=1:16
    if  isequal(B(i:i+1,j),C(:,1,k))==1  %isequal(C(:,1,RandV(temp)),C(:,1,k))==1 
n=n+1;
k_col(n)=k;
    end
end

freq=count(k_col); 

%----Pick Random configurations based on the subset of the Histogram -----
Randfreq=RandHistValues(freq,k_col);
 
 
    B(i:i+1,j+1)=C(:,2,Randfreq(temp));
    
end
end


%---REST OF THE NEW IMAGE-----
clear freq
clear Randfreq


for i=2:c:length(B(:,1))-1   
for j=1:c:length(B)-1
        

 
   if j==1 %mod(j,2)==1 %If j is odd
 
        %tempodd=tempodd+1;
  %--STEP2: Find patterns whose 1st ROW is the same as the current pattern---
n=0;
for k=1:16
    if  isequal(B(i,j:j+1),C(1,:,k))==1 
n=n+1;
k_row(n)=k;
    end
end

freq=count(k_row); 

%----Pick Random configurations based on the subset of the Histogram -----
Randfreq=RandHistValues(freq,k_col);
 

    B(i+1,j:j+1)=C(2,:,Randfreq(1));
 
    
    
    
   else  %j>1  %j is even
       %tempeven=tempeven+1;
    %--STEP3: Find 2x2 patterns whose unique unknown pixel is the botton rigth one---   
       n=0;
for k=1:16
    if  isequal(B(i,j:j+1),C(1,:,k))==1  &&  isequal(B(i+1,j),C(2,1,k))==1
n=n+1;
k_3(n)=k;
    end
end

freq3=count(k_3); 
     
%----Pick Random configurations based on the subset of the Histogram -----
Randfreq3=RandHistValues(freq3,k_3);
 

    B(i+1,j+1)=C(2,2,Randfreq3(1));

    
   end 
   
  
    
end
end





figure(4)
%image(B)    % Here white looks yellow and black looks blue (if not binary)
imshow(B)    %Show the black- white image





