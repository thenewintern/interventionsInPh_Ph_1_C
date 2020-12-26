%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PhPh1c_PE_method_2 file. Takes in the current 0th, 1st and 2nd conditional
%moments along with the change (magnitude) and choice (addition or
%deletion). This function calculates the required probability values by
%successive conditioning. The logic is fairly simple. Since the PE
%approximation works well while approximating boundary probabilities, we
%try to approximate each probability value by treating it as a boundary
%probability. Let us illustrate this with an example. Lets say c=10 and we
%wish to add 3 entities to the system. Thus we would need the probabilities
%for being in state 8,9 and 10 for all possible arrival and service phase
%combinations. We first approximate the probability of being in state 10.
%Once we have this, we subtract this value from the 0th conditional moment,
%c*value from the first conditional moment and c^2*value from the second 
%conditional moment. Thus for the current set of values of the conditional
%moments, state 9 is a boundary probability. We extract this value and
%perform the procedure again to ensure the next time we can extract the
%probability of being in state 8 as a boundary probability. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = PhPh1c_PE_method_2(E0,E1,E2,change,choice)


global c

P=zeros(c,1);
den=c;


if choice==1

count=c;    
    
for k=c-change+1:c

    if E0>0&&E1>0&&E2>0
        
            EX=E1/E0;
            EX2=E2/E0;
    
            EN=EX-1;
            EN2=EX2-2*EN-1;
        
    
    if isnan(EN)~=1 
           
        den=den-1;
           
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

          
          P(count,1)=PhPh1c_PE(den,1-theta,gamma);
          if P(count,1)<0
              P(count,1)=0;
          end
          
          
    end
    
    P(count,1)=P(count,1)*E0;
    %Update values of the conditional moments
    E0=E0-P(count,1);
    E1=E1-count*P(count,1);
    E2=E2-(count^2)*P(count,1);
    count=count-1;
    
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if choice==2
    
    
for k=1:change

    if E0>0&&E1>0&&E2>0
        
            EX=E1/E0;
            EX2=E2/E0;
    
            EN=EX-k;
            EN2=EX2-2*k*EN-(k^2);
        
    
    if isnan(EN)~=1 
           
        den=den-1;
           
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

          
          P(k,1)=PhPh1c_PE(den,theta,gamma);
          if P(k,1)<0
              P(k,1)=0;
          end
          
          
    end
    
    P(k,1)=P(k,1)*E0;
    %Update values of the conditional moments
    E0=E0-P(k,1);
    E1=E1-k*P(k,1);
    E2=E2-(k^2)*P(k,1);
    
    end
end

end


end
