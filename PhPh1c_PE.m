%Polya Eggenberger approximation for the boundary probabilities. 
function [x]=PhPh1c_PE(den,theta,gamma)
    frac = 1;
    for m = 1:den
        n=m-1;
        frac = frac*(1-theta+n*gamma)/(1+n*gamma);
    end
    x=frac;
end
    