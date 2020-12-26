%%%%%%%%%% Arrival and Service related matrices and vectors %%%%%%%%%%%%%%%
% function[k1,k2,alpha,beta,A1,A2,B1,B2,lambda,mu]=PhPh1c_qparm(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k1,k2,alpha,beta,A1,A2,B1,B2,lambda,mu]=PhPh1c_qparm(t)

%All parameters are currently constants. The user may make the parameters
%dependent on time by using t (That is the input argument for this function
%when it is called by PhPh1c_function file). 

%Vector names and usage :
%lambda : row vector consisting of arrival rates for each arrival phase
%mu : row vector consisting of service rates for each service phase
%alpha : row vector consisting of initial arrival phase probabilities
%beta : row vector consisting of initial service phase probabilities
%A1 : square matrix consisting of transition probabilities for arrival
%phases
%B1 : square matrix consisting of transition probabilities for service 
%phases
%A2 : arrival phase absorption probabilities 
%B2 : service phase absorption probabilities

%All matrices are set for 10 phases. In order to change this, just change
%the size of the matrices to the maximum possible number of phases that can
%be associated with them. 

lambda=[4  0  0  0  0  0  0  0  0  0];

mu=[5  0  0  0  0  0  0  0  0  0];


alpha=[1.  0.  0.  0.  0.  0.  0.  0  0  0];

beta=[1.  0.  0  0  0  0  0  0  0  0];


A1=[0.  0.  0.  0  0.  0  0  0  0  0;   %1
    
    0.  0.  0.  0  0  0  0  0  0  0;   %2
    
    0.  0.  0.  0  0  0  0  0  0  0;   %3
    
    0  0  0  0  0  0  0  0  0  0;   %4
    
    0  0  0  0  0  0  0  0  0  0;   %5
    
    0  0  0  0  0  0  0  0  0  0;   %6
    
    0  0  0  0  0  0  0  0  0  0;   %7
    
    0  0  0  0  0  0  0  0  0  0;   %8
    
    0  0  0  0  0  0  0  0  0  0;   %9
    
    0  0  0  0  0  0  0  0  0  0];  %10


B1=[0.  0.  0  0  0  0  0  0  0  0;   %1
    
    0.  0.  0  0  0  0  0  0  0  0;   %2
    
    0  0  0  0  0  0  0  0  0  0;   %3
    
    0  0  0  0  0  0  0  0  0  0;   %4
    
    0  0  0  0  0  0  0  0  0  0;   %5
    
    0  0  0  0  0  0  0  0  0  0;   %6
    
    0  0  0  0  0  0  0  0  0  0;   %7
    
    0  0  0  0  0  0  0  0  0  0;   %8
    
    0  0  0  0  0  0  0  0  0  0;   %9
    
    0  0  0  0  0  0  0  0  0  0];  %10


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section eliminates as many 0 entries in the matrices as possible. 
%This is done by checking for the non zero element with the highest index 
%in each of the vectors lambda, mu, alpha and beta along with the matrices 
%A1 and B1. The logic is, beyond this index no value contributes to a 
%change in the differential equations. The highest index for arrival 
%processes related vectors (lambda and alpha) and matrix (A1) that 
%corresponds to a non zero element is then designated as the number of 
%arrival phases k1. A similar process with mu, beta and B1 determines 
%the number of service phases k2. To see the implementation of this logic,
%refer to the PhPh1c_function file. 
sum_alpha=sum(alpha);
sum_beta=sum(beta);

if sum_alpha~=1
    alpha=alpha/sum_alpha;
end

if sum_beta~=1
    beta=beta/sum_beta;
end

%k1=min(max(find(lambda)),max(find(alpha)));
%k2=min(max(find(mu)),max(find(beta)));

k1=max(find(lambda));
k2=max(find(mu));

A2=1-(sum(A1'))';
B2=1-(sum(B1'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

