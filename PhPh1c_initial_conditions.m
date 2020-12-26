%PhPh1c_initial_conditions function file. Has three input arguments. The
%first argument is to choose whether to make a deterministic change to the
%system (choice=1) or not (choice=0). If choice=0, then the function has
%been invoked at the start of the run when the first set of initial
%conditions are needed. Thereafter, everytime the function is called, it is
%to make a deterministic change.
%The second argument is the deterministic change to be made. This value is
%taken from the det_changes vector.
%The third argument is a vector of size g1+c*g1*g2+3*g1*g2+g1. At the end
%of every run in the PhPh1c_det_changes file (Each run is from time of
%previous change to the time of next change in the system), the last row of
%values obtained by running the ode45 is taken as the initial conditions
%for the next run. Since we have to make deterministic changes to the
%system at that time, these values must be changed accordingly. The third
%argument for this function is the last set of values obtained from the
%previous run. This is used as a column vector x. Note, the size of this
%vector is g1+c*g1*g2+3*g1*g2+g1 for the first g1+c*g1*g2 values are the
%values of the state probabilities. The next 3*g1*g2 values are for the 3
%sets of PMDEs and the last g1 values are for the g1 idle states. Note: 
%g1 = max_arrival_phases and g2 = max_service_phases. 
%The program outputs the vector named ics (initial conditions) that is used
%as the initial conditions for the next run. 

%%%%%%%%%%%%%%%%%%%%%% Initial conditions file %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ics]=PhPh1c_initial_conditions(choice,change,x)

global c max_arrival_phases max_service_phases;

epsilon=1e-14;

ics=zeros(max_arrival_phases+...
           (c+3)*max_arrival_phases*max_service_phases+...
            max_arrival_phases,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%When choice = 0, the logic is fairly simple and thus will not be explained
%here. 

if choice==0


initial_number_of_entities=change;

if initial_number_of_entities==0
    [k1,k2,alpha,]=PhPh1c_qparm(0);
    
    ics(1:max_arrival_phases,1)=alpha(1,1:max_arrival_phases)';
    ics(max_arrival_phases+...
        (c+3)*max_arrival_phases*max_service_phases+1:...
        max_arrival_phases+...
        (c+3)*max_arrival_phases*max_service_phases+max_arrival_phases...
        ,1)=alpha(1,1:max_arrival_phases)';
else

if initial_number_of_entities>0&&initial_number_of_entities<=c

    [k1,k2,alpha,beta,]=PhPh1c_qparm(0);
    
    dummy=alpha'*beta;
    
    for i=1:max_arrival_phases
        for j=1:max_service_phases
            index=max_arrival_phases+(initial_number_of_entities-1)*...
                                  max_arrival_phases*max_service_phases+...
                                  (i-1)*max_service_phases+j;
            
            E0_index=max_arrival_phases+c*max_arrival_phases*...
                     max_service_phases+(i-1)*max_service_phases+j;
                 
            E1_index=max_arrival_phases+(c+1)*max_arrival_phases*...
                     max_service_phases+(i-1)*max_service_phases+j;
            
            E2_index=max_arrival_phases+(c+2)*max_arrival_phases*...
                     max_service_phases+(i-1)*max_service_phases+j;
                 
                 
            ics(index)=dummy(i,j);
            ics(E0_index)=dummy(i,j);
            ics(E1_index)=initial_number_of_entities*dummy(i,j);
            ics(E2_index)=(initial_number_of_entities^2)*dummy(i,j);
            
        end
    end
       
else
    [k1,k2,alpha,]=PhPh1c_qparm(0);
    
   ics(1:max_arrival_phases,1)=alpha(1,1:max_arrival_phases)';
    ics(max_arrival_phases+...
        (c+3)*max_arrival_phases*max_service_phases+1:...
        max_arrival_phases+...
        (c+3)*max_arrival_phases*max_service_phases+max_arrival_phases...
        ,1)=alpha(1,1:max_arrival_phases)';
end

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The heart of the code lies in the next section which handles deterministic
%changes. We can have two types of deterministic changes, positive or
%negative. Both will be explained.


%First we handle the Kolmogorov Forward equations. The logic for which is
%fairly simple. For deterministic additions, we map lower states to higher
%states. The higher state is determined by the magnitude of the addition.
%Note that we can map a lower state to a higher state for 0,1,...,c-change,
%where change denotes the deterministic addition. For states higher than
%c-change, we map it to state c. Mapping essentially means shifting the
%current probability values from lower states to higher states. And for
%states greater than c-change, we 'accumulate' it in state c. The logic is
%the same for deterministic addition, except now we map higher states to
%lower states. For the PMDEs, let us explain it with an example of
%deterministic addition. All probability values for states 0,1,...,c-change
%get shifted up by the magnitude of the change. But it's the states
%c-change+1 to c that matter to us, since they get mapped to state c.
%Hence, we would need the current probability values for states
%c-change+1,...,c. Instead of solving the Kolmogorov Forward equations for
%this, we 'guess' these probabilities (The programmer wishes to point out
%the irony in guessing probability values, when probability itself can be
%described as a guess. In this regard, the programmer wishes to apologize
%for the poor choice of words) using the Polya Eggenberger distribution.
%There are two methods to do this, Method 1 involves using the distribution
%as it is and obtaining the different probability values. Method 2 relies
%on the uncanny nature of the PE distribution to produce good guesses of
%boundary probabilities. So in Method 2, we successively condition and
%obtain the required probabilities as boundary probabilities. To know more
%about this, refer to PhPh1c_PE_method_1 and PhPh1c_PE_method_2 files. 
if choice==1

    [k1,k2,alpha,beta,]=PhPh1c_qparm(0);  
    
    %declare all variables.
    
    E0_pmde=zeros(max_arrival_phases,max_service_phases);
    E1_pmde=zeros(max_arrival_phases,max_service_phases);
    E2_pmde=zeros(max_arrival_phases,max_service_phases);
    
    E0_dummy=zeros(max_arrival_phases,max_service_phases);
    E1_dummy=zeros(max_arrival_phases,max_service_phases);
    E2_dummy=zeros(max_arrival_phases,max_service_phases);
    
    
    P0_feq=x(1:max_arrival_phases,1); %KFE values for the idle states
    
    P0_PE=x(max_arrival_phases+...
     (c+3)*max_arrival_phases*max_service_phases+1:...
     max_arrival_phases+...
     (c+3)*max_arrival_phases*max_service_phases+max_arrival_phases,1); %Idle state values obtained by using the PE approximation

    P_feq=zeros(max_arrival_phases,max_service_phases,c);
    P_dummy_feq=zeros(max_arrival_phases,max_service_phases,c);
       
    
for n=1:c
    for i=1:max_arrival_phases
        for j=1:max_service_phases
            index=max_arrival_phases+...
                  (n-1)*max_arrival_phases*max_service_phases+...
                  (i-1)*max_service_phases+j;
              
             P_dummy_feq(i,j,n)=x(index,1); %All non idle state probabilities
             %are stored in a dummy vector. The use of the dummy vector
             %will soon become evident. 
              
        end
    end
end
   
for i=1:max_arrival_phases
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
             
        
        %All current PMDE values are stored in dummy vectors.      
        E0_dummy(i,j)=x(E0_index,1); 
        E1_dummy(i,j)=x(E1_index,1);
        E2_dummy(i,j)=x(E2_index,1);
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if change>0
    
    if change>c
        change=c;
    end
    
%State 0 gets mapped to the state with the value of |change|
P0_feq=P0_feq*beta;

P_feq(:,:,change)=P0_feq(:,:);

%Map all states accordingly
for n=1:c-change
    P_feq(:,:,n+change)=P_dummy_feq(:,:,n);
end

for n=c-change+1:c
    P_feq(:,:,c)=P_feq(:,:,c)+P_dummy_feq(:,:,n);
end

%Once the states are mapped, store them in the ics vector. 
for n=1:c
    for i=1:max_arrival_phases
        for j=1:max_service_phases
            index=max_arrival_phases+...
                  (n-1)*max_arrival_phases*max_service_phases+...
                  (i-1)*max_service_phases+j;
              
            ics(index,1)=P_feq(i,j,n);
            
        end
    end
end

%Do the same with the idle state probability values calculated using the PE
%approximation. 
P0_PE=P0_PE*beta;

%We need the probability values for the states c-change+1 to c, and all we
%have as information are the 0th, 1st and 2nd PMDE. Use either method 1 or
%method 2 to obtain the required probability values. 
for i=1:max_arrival_phases
    for j=1:max_service_phases
        
        
        P_dummy_PE=...
        PhPh1c_PE_method_2(E0_dummy(i,j),E1_dummy(i,j),E2_dummy(i,j),...
                           change,1);
        
        %Once we've obtained the probability values, begin updating all 
        %PMDE values for the next run. 
        E0_pmde(i,j)=E0_dummy(i,j)+P0_PE(i,j);
        
        PN=0;
        PN2=0;
        sum_PN=0;
        dummy1=0;
        dummy2=0;
        
        multiplier=change-1;
        
        for k=c-change+1:c
            sum_PN=sum_PN+P_dummy_PE(k,1);
            
            PN=PN+P_dummy_PE(k,1)*multiplier;
            PN2=PN2+P_dummy_PE(k,1)*(multiplier^2);
            
            dummy1=dummy1+k*P_dummy_PE(k,1);
            dummy2=dummy2+k*multiplier*P_dummy_PE(k,1);
            
            multiplier=multiplier-1;
        end
        
        E1_pmde(i,j)=E1_dummy(i,j)+change*P0_PE(i,j)+...
                (E0_dummy(i,j)-sum_PN)*change+PN;        
        
        E2_pmde(i,j)=E2_dummy(i,j)+(change^2)*P0_PE(i,j)+...
                2*((E1_dummy(i,j)-dummy1)*change+dummy2)+...
                (E0_dummy(i,j)-sum_PN)*(change^2)+PN2;
            
    end
end

%Once all the calculations are performed, store all values in the ics
%vector accordingly. 
for i=1:max_arrival_phases
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
        
        
        ics(E0_index,1)=E0_pmde(i,j);     
        ics(E1_index,1)=E1_pmde(i,j);
        ics(E2_index,1)=E2_pmde(i,j);
        
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The logic is similar when the change made to the system is negative. In
%this case it becomes tricky since we have to map probability values to the
%idle state which has only arrival phases. Thus for higher states we sum up
%all probability values for one particular arrival phases in all possible 
%service phases and map it to the same arrival phase in the idle state. 

if change<0
    
  change=abs(change);
  
  if change>c
      change=c;
  end
  
  %Map higher states to lower states for the KFE
  for n=1:change
      for i=1:max_arrival_phases
          for j=1:max_service_phases
              P0_feq(i,1)=P0_feq(i,1)+P_dummy_feq(i,j,n);              
          end
      end
  end
  
  for n=change+1:c
      for i=1:max_arrival_phases
          for j=1:max_service_phases
              P_feq(i,j,n-change)=P_dummy_feq(i,j,n);
          end
      end
  end
  
  ics(1:max_arrival_phases,1)=P0_feq(:,1);
  
  for n=1:c
      for i=1:max_arrival_phases
          for j=1:max_service_phases
              index=max_arrival_phases+...
                    (n-1)*max_arrival_phases*max_service_phases+...
                    (i-1)*max_service_phases+j;
                
              ics(index,1)=P_feq(i,j,n);
              
          end
      end
  end  
  
  %Once we're done with the KFE, we now shift our attention to the second
  %set of equations. In this case, we will have to guess the probability
  %values of the required states and add these to the corresponding idle
  %state.
  
  for i=1:max_arrival_phases
      for j=1:max_service_phases
        
        P_dummy_PE=...
        PhPh1c_PE_method_2(E0_dummy(i,j),E1_dummy(i,j),E2_dummy(i,j),...
                           change,2);
          
       
        sum_PN=0;
        PN=0;
        PN2=0;
        
        for n=1:change
            sum_PN=sum_PN+P_dummy_PE(n,1);
            PN=PN+n*P_dummy_PE(n,1);
            PN2=PN2+(n^2)*P_dummy_PE(n,1);
            P0_PE(i,1)=P0_PE(i,1)+P_dummy_PE(n,1);
        end
        
        E0_pmde(i,j)=E0_dummy(i,j)-sum_PN;
        
        E1_remaining=E1_dummy(i,j)-PN;
        E2_remaining=E2_dummy(i,j)-PN2;
        
        E1_pmde(i,j)=E1_remaining-change*E0_pmde(i,j);
        E2_pmde(i,j)=E2_remaining+(change^2)*E0_pmde(i,j)-...
                     2*change*E1_remaining;
                 
        if E0_pmde(i,j)<epsilon
            E0_pmde(i,j)=0;
        end
        
        if E1_pmde(i,j)<epsilon
            E1_pmde(i,j)=0;
        end
        
        if E2_pmde(i,j)<epsilon
            E2_pmde(i,j)=0;
        end
        
        
      end
  end
 
  %Once we're done mapping all states, store in the ics array. 
  
  for i=1:max_arrival_phases
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
        
        ics(E0_index,1)=E0_pmde(i,j);     
        ics(E1_index,1)=E1_pmde(i,j);     
        ics(E2_index,1)=E2_pmde(i,j);     
      end
  end
  
  ics(max_arrival_phases+(c+3)*max_arrival_phases*max_service_phases+1:...
      max_arrival_phases+(c+3)*max_arrival_phases*max_service_phases+...
      max_arrival_phases,1)=P0_PE(:,1);
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end