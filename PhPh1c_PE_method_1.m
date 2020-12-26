%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PhPh1c_PE_method_1 file. Takes in the current 0th, 1st and 2nd conditional
%moments and then returns the probability values for states 0,1,...,c-1.
%Note the input arguments contain two dummy variables. These are not used
%by the code. They have been left there to allow the user to choose any
%method in the PhPh1c_initial_conditions file by just changing the number
%(1 or 2) in the function name. PhPh1c_PE_method_2 function does use 5
%arguments. For more information, refer to the PhPh1c_PE_method_2 function
%file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P] = PhPh1c_PE_method_1(E0,E1,E2,dummy1,dummy2)


global c
c=10;

P=zeros(c,1);

            EX=E1/E0;
            EX2=E2/E0;
    
            EN=EX-1;
            EN2=EX2-2*EN-1;
        
    
    if isnan(EN)~=1 
           
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

          
     
frac = 1.;
for j = 0:den-1;
    frac = frac * (1. - theta + j * gamma) / (1. + j * gamma);
 end
P(1,1) =  frac;
%
% PX(C+1) = P(X = C)
%
frac = 1.;
for j = 0:den-1;
    frac = frac * ( theta + j * gamma) / (1. + j * gamma);
 end
P(c,1) =  frac;
%
% PX(i+1) = P(X = i), for i = 1,c-1
%
for i = 1:den-1
    frac = 1.;
    for j = 0:i-1;
        frac = frac * (theta + j * gamma)/(1. + j * gamma);
    end
    for j = i:c-1;
        frac = frac * ( (1. - theta) + (j - i) * gamma ) / ( 1. + j * gamma );
    end
    P(i+1) = nchoosek(den,i) * frac;
end
    
    end
    
    P=P*E0;
    
end


