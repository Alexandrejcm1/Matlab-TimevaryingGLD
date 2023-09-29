function [r,e2,e3,e4,e5]=simulTukey_lambda3t(draws,sampleSize,burn,lambda1,lambda4,theta,omega2,alpha2,beta2,omega3,alpha3,beta3);
draws = 1;
sampleSize=550;
burn  = 50;

% imposing lambda1=0, lambda3=lambda4 => please generalize for any lambda

alpha2 = 0.05;
beta2  = 0.9;
alpha3 = 0.01;
beta3  = 0.9;
omega2 = (1/(2*lambda3+1)- 1/(2*lambda3^2+2*lambda3+1)+2/(lambda3*lambda4+lambda3+lambda4+1)-1/(2*lambda4^2+2*lambda4+1)+1/(2*lambda4+1)-beta(lambda3+1,lambda4+1))*(1-alpha2-beta2)/sigma2;
omega3 = (1/(((1/2*lambda3+1)+(1/2*lambda4+1)-2*beta(lambda3+1,lambda4+1))-((1/lambda3+1)-(1/lambda4+1))^2)^1.5*((1/3*lambda3+1)-(1/3*lambda4+1)-3*beta(2*lambda3+1,lambda4+1)+3*beta(lambda3+1,2*lambda4+1)))-(3*((1/(lambda3+1)-(1/lambda4+1))*((1/(2*lambda3+1)+(1/2*lambda4+1)-2*beta(lambda3+1,lambda4+1)))))+(2*((1/lambda3+1)-(1/lambda4+1))^3)*(1-alpha3-beta3)/sigma3;

theta = 0.95;

lambda2 = omega2/(1-alpha2-beta2);
u0 = rand(1,draws);
r0 = lambda1+((u0.^lambda3-(1-u0).^lambda4)/(-sqrt(lambda2)));
u1 = rand(1,draws);

lambda2t = ones(1,draws).*(omega2+alpha2*r0.^2+beta2*lambda2);
ret      = ones(1,draws).*(lambda1+(u1.^lambda3-(1-u1).^lambda4)/(-sqrt(omega2+alpha2*r0.^2+beta2*lambda2t)));
ewma2    = lambda2t;
ewma3    = zeros(1,draws);
ewma4    = zeros(1,draws);
ewma5    = zeros(1,draws);
for j=1:draws
   for i=1:sampleSize-1
       lambda2t(i+1,j) = omega2+alpha2*ret(i,j)^2+beta2*lambda2t(i,j);
       u               = rand(1,1);
       ret(i+1,j)      = lambda1+(u.^lambda3-(1-u).^lambda4)/(-sqrt(lambda2t(i+1,j)));
       ewma2(i+1,j)    = (1-theta)*ret(i,j)^2+theta*ewma2(i,j);
       ewma3(i+1,j)    = (1-theta)*ret(i,j)^3+theta*ewma3(i,j);
       ewma4(i+1,j)    = (1-theta)*ret(i,j)^4+theta*ewma4(i,j);
       ewma5(i+1,j)    = (1-theta)*ret(i,j)^5+theta*ewma5(i,j);
   end
   % burning initial 50 observations
   r(:,j)   = ret(51:sampleSize,j);
   e2(:,j) = ewma2(51:sampleSize,j);
   e3(:,j) = ewma3(51:sampleSize,j);
   e4(:,j) = ewma4(51:sampleSize,j);
   e5(:,j) = ewma5(51:sampleSize,j);
end

  
  
end



