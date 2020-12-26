%%%%%%%%%%%%%%%%%%% Differential equations file %%%%%%%%%%%%%%%%%%%%%%%%%%%
%PhPh1c_function file. This file contains the implementation of the
%differential equations (Kolmogorv Forward Equations and the PMDE). For
%detailed information on the equations implemented, refer to the paper :
%Approximating Nonstationary Ph(t)/Ph(t)/1/c Queueing Systems (Ong and
%Taaffe - 1988)
%For more information on the logic behind the differential equations and
%the Polya Eggenberger approximation technique, refer to
%the research paper mentioned above. 
function [dpdt]=PhPh1c_function(t,y)
t %This is ornamental and was left here for the programmer to keep track of
%the calculations. This may be removed by the user. 
global c max_arrival_phases max_service_phases;

%The ode45 will return two arrays. t contains the values of time for which
%the ode45 calculates values for all the differential equations set in this
%file. y is returned as x (Refer to PhPh1c_det_changes to check this fact).
%y essentially consists of two sets of differential equations. The first
%set is the Kolmogorov forward equations.
%max_arrival_phases+c*max_arrival_phases*max_service_phases denotes the
%number of Kolmogorov Forward equations that are calculated. The next
%3*max_arrival_phases*max_service_phases+max_arrival_phases are arranged in
%the following way: The first max_arrival_phases*max_service_phases number
%of equations are the 0th PMDEs for the Ph(t)/Ph(t)/1/c queue, the next
%max_arrival_phases*max_service_phases number of equations are the 1st PMDE
%equations and the following max_arrival_phases*max_service_phases number
%of equations are the 2nd PMDEs. Since we need the probability values for
%each arrival phase when the system is idle, the last max_arrival_phases
%number of equations represent the differential equations of the
%Ph(t)/Ph(t)/1/c system when its idle. In order to calculate the PMDEs and
%the idle state differential equations we need the values for P_{1,i,j} 
%and P_{c,i,j}. This is accomplished by using the PE approximation.  

dpdt=zeros(max_arrival_phases+...
           (c+3)*max_arrival_phases*max_service_phases+...
            max_arrival_phases,1);

[k1,k2,alpha,beta,A1,A2,B1,B2,lambda,mu]=PhPh1c_qparm(t);



P=zeros(max_arrival_phases,max_service_phases,c);  %stores all current probability values of the non idle states
E0=zeros(max_arrival_phases,max_service_phases,c); %stores the current values of the 0th moment PMDE values
E1=zeros(max_arrival_phases,max_service_phases,c); %stores the current values of the 1st moment PMDE values
E2=zeros(max_arrival_phases,max_service_phases,c); %stores the current values of the 2nd moment PMDE values
P1=zeros(max_arrival_phases,max_service_phases);   %stores the approximation values of all possible arrival and service phase combination when there is one entity in the queue
PC=zeros(max_arrival_phases,max_service_phases);   %stores the approximation values of all possible arrival and service phase combination when there are c entities in the queue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P0=y(1:max_arrival_phases,1);

for n=1:c
    for i=1:max_arrival_phases
        for j=1:max_service_phases
            P(i,j,n)=y(max_arrival_phases+(n-1)*max_arrival_phases*...
                                                max_service_phases+...
                                                (i-1)*max_service_phases...
                                                +j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Forward Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Refer to the research paper mentioned above to understand the logic and
%the terminology associated wih the following equations. 

for i=1:max_arrival_phases
    index=i;
    
    dpdt(index)=-lambda(1,i)*P0(i,1)+sum(A1(1:k1,i)'.*lambda(1,1:k1).*...
                                          P0(1:k1,1)')+...
                         sum(B2(1:k2,1)'.*mu(1,1:k2).*P(i,1:k2,1));
                     
end

for n=1:c
    for i=1:max_arrival_phases
        for j=1:max_service_phases
            index=max_arrival_phases+(n-1)*max_arrival_phases*...
                                     max_service_phases+...
                                     (i-1)*max_service_phases+j;
                                 
            
            dpdt(index)=-(lambda(1,i)+mu(1,j))*P(i,j,n)+...
                         sum(A1(1:k1,i)'.*lambda(1,1:k1).*P(1:k1,j,n)')+...
                         sum(B1(1:k2,j)'.*mu(1,1:k2).*P(i,1:k2,n));
                     
            
            if n==1
                dpdt(index)=dpdt(index)+...
                            alpha(1,i)*beta(1,j)*...
                            sum(A2(1:k1,1)'.*lambda(1,1:k1).*P0(1:k1,1)');
            end
            
            if n>1
                dpdt(index)=dpdt(index)+alpha(1,i)*...
                          sum(A2(1:k1,1)'.*lambda(1,1:k1).*P(1:k1,j,n-1)');
                      
                      if n==c
                          dpdt(index)=dpdt(index)+alpha(1,i)*...
                                      sum(A2(1:k1,1)'.*lambda(1,1:k1).*...
                                          P(1:k1,j,c)');
                      end
                      
            end
            
            if n<c
                dpdt(index)=dpdt(index)+beta(1,j)*...
                          sum(B2(1:k2,1)'.*mu(1,1:k2).*P(i,1:k2,n+1));
            end
            
        end
    end
       
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear P0 P;

%%%%%%%%%%%%%%%%%%%%% Polya Approximation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max_arrival_phases
    
    P0(1,i)=y(max_arrival_phases+...
              (c+3)*max_arrival_phases*max_service_phases+i,1);
              
    
    for j=1:max_service_phases
        
        E0(i,j)=y(max_arrival_phases+...
                  c*max_arrival_phases*max_service_phases+...
                  (i-1)*max_service_phases+j,1);
        
        E1(i,j)=y(max_arrival_phases+...
                  (c+1)*max_arrival_phases*max_service_phases+...
                  (i-1)*max_service_phases+j,1);
        
        E2(i,j)=y(max_arrival_phases+...
                  (c+2)*max_arrival_phases*max_service_phases+...
                  (i-1)*max_service_phases+j,1);
              
        
              
            EX=E1(i,j)/E0(i,j);
            EX2=E2(i,j)/E0(i,j);
    
            %%%%%%%%%%%%%%%%%%%%%
            EN=EX-1;
            EN2=EX2-2*EN-1;
            %%%%%%%%%%%%%%%%%%%%%%
            %The previous two lines of code are not immediately obvious.
            %Note that when we approximate P_{1,i,j} and P_{c,i,j} using
            %E0(i,j),E1(i,j) and E2(i,j) we do so with a state space that
            %consistes of indices 0,1,...,c-1. So we have to shift our
            %conditional moments appropriately to obtain P_{1,i,j} and
            %P_{c,i,j}
        
    
    if isnan(EN)==1 
        P1(i,j)=0;
        PC(i,j)=0;
        E1(i,j)=0;
        E2(i,j)=0;
        E0(i,j)=0;
    else
        
        % This section of the code has been provided by Dr. Michael
        % Taaffe,ISE Department, Virignia Tech.
        %%%%%     %%%%%     %%%%%     %%%%%     %%%%%     %%%%%     %%%%%
        den=c-1;
           
        Var=EN2-EN*EN;
        eita=1000.;
        x=EN/den;
        epsilon=10^(-4);
          if EN<epsilon
              p=0.;
          else
              p=x;
          end
        theta=p;

        q=1-p;
        pqmin=min(p,q);
        x=-pqmin/(den-1);
        d=EN2-den*EN;
          if(EN<epsilon) 
              gamma=0.;
          elseif(EN>(den-epsilon))
              gamma=0.;
          elseif(Var<epsilon)
              gamma=-1/den+epsilon;
          elseif(abs(d)<epsilon)
              gamma=eita;
          elseif ((EN*(EN+(1.-p)))-EN2)/d<x
              gamma=x+epsilon; 
          else
              gamma=((EN*(EN+(1.-p)))-EN2)/d;
          end
          
          %%%%%     %%%%%     %%%%%     %%%%%     %%%%%     %%%%%     %%%%%

          P1(i,j)=PhPh1c_PE(den,theta,gamma);
          PC(i,j)=PhPh1c_PE(den,1-theta,gamma);
          
          P1(i,j)=P1(i,j)*E0(i,j);
          PC(i,j)=PC(i,j)*E0(i,j);
          
    end
    
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%% Second Set of Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max_arrival_phases
    
    P0_index=max_arrival_phases+...
              (c+3)*max_arrival_phases*max_service_phases+i;
           
    dpdt(P0_index)=-lambda(1,i)*P0(1,i)+...
                   sum(A1(1:k1,i)'.*lambda(1,1:k1).*P0(1,1:k1))+...
                   sum(B2(1:k2,1)'.*mu(1,1:k2).*P1(i,1:k2));
    
    
    for j=1:max_service_phases
        
          
        E0_index=max_arrival_phases+...
                 c*max_arrival_phases*max_service_phases+...
                 (i-1)*max_service_phases+j;
             
        E1_index=max_arrival_phases+...
                 (c+1)*max_arrival_phases*max_service_phases+...
                 (i-1)*max_service_phases+j;
             
        E2_index=max_arrival_phases+...
                 (c+2)*max_arrival_phases*max_service_phases+...
                 (i-1)*max_service_phases+j;
        
                     
        dpdt(E0_index)=-(lambda(1,i)+mu(1,j))*E0(i,j)+...
                       sum(A1(1:k1,i)'.*lambda(1,1:k1).*E0(1:k1,j)')+...
                       sum(B1(1:k2,j)'.*mu(1,1:k2).*E0(i,1:k2))+...
                       beta(1,j)*(sum(B2(1:k2,1)'.*mu(1,1:k2).*E0(i,1:k2))...
                       -sum(B2(1:k2,1)'.*mu(1,1:k2).*P1(i,1:k2)))+...
                       alpha(1,i)*sum(A2(1:k1,1)'.*lambda(1,1:k1).*E0(1:k1,j)')+...
                       alpha(1,i)*beta(1,j)*...
                       sum(A2(1:k1,1)'.*lambda(1,1:k1).*P0(1,1:k1));
        
        dpdt(E1_index)=-(lambda(1,i)+mu(1,j))*E1(i,j)+...
                       sum(A1(1:k1,i)'.*lambda(1,1:k1).*E1(1:k1,j)')+...
                       sum(B1(1:k2,j)'.*mu(1,1:k2).*E1(i,1:k2))+...
                       beta(1,j)*(sum(B2(1:k2,1)'.*mu(1,1:k2).*E1(i,1:k2))-...
                       sum(B2(1:k2,1)'.*mu(1,1:k2).*E0(i,1:k2)))+...
                       alpha(1,i)*(sum(A2(1:k1,1)'.*lambda(1,1:k1).*E1(1:k1,j)')+...
                       sum(A2(1:k1,1)'.*lambda(1,1:k1).*E0(1:k1,j)')-...
                       sum(A2(1:k1,1)'.*lambda(1,1:k1).*PC(1:k1,j)'))+...
                       alpha(1,i)*beta(1,j)*sum(A2(1:k1,1)'.*lambda(1,1:k1).*P0(1,1:k1));
                   
        dpdt(E2_index)=-(lambda(1,i)+mu(1,j))*E2(i,j)+...
                       sum(A1(1:k1,i)'.*lambda(1,1:k1).*E2(1:k1,j)')+...
                       sum(B1(1:k2,j)'.*mu(1,1:k2).*E2(i,1:k2))+...
                       beta(1,j)*(sum(B2(1:k2,1)'.*mu(1,1:k2).*E2(i,1:k2))-...
                       2*sum(B2(1:k2,1)'.*mu(1,1:k2).*E1(i,1:k2))+...
                       sum(B2(1:k2,1)'.*mu(1,1:k2).*E0(i,1:k2)))+...
                       alpha(1,i)*(sum(A2(1:k1,1)'.*lambda(1,1:k1).*E2(1:k1,j)')+...
                       2*sum(A2(1:k1,1)'.*lambda(1,1:k1).*E1(1:k1,j)')+...
                       sum(A2(1:k1,1)'.*lambda(1,1:k1).*E0(1:k1,j)')-...
                       (2*c+1)*sum(A2(1:k1,1)'.*lambda(1,1:k1).*PC(1:k1,j)'))+...
                       alpha(1,i)*beta(1,j)*sum(A2(1:k1,1)'.*lambda(1,1:k1).*P0(1,1:k1));
                   
                   
    end
    
end
                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end