function [r,e2,e3,e4]=simulTukey(sampleSize,burn,lambda4,theta,lambda3,sigma2,alpha2,beta2);
draws = 2000;

A = 1/(lambda3+1)-1/(lambda4+1);
B = 1/(2*lambda3+1)+1/(2*lambda4+1)-2*beta(lambda3+1,lambda4+1);
C = 1/(3*lambda3+1)-1/(3*lambda4+1)-3*beta(2*lambda3+1,lambda4+1)+3*beta(lambda3+1,2*lambda4+1);
D = 1/(4*lambda3+1)+1/(4*lambda4+1)-4*beta(3*lambda3+1,lambda4+1)+6*beta(2*lambda3+1,2*lambda4+1)-4*beta(lambda3+1,3*lambda4+1);

lambda2 = sqrt((B-A^2)/sigma2);
omega2  = (1-beta2)*lambda2-alpha2*(B-A^2);
lambda1 = -A/lambda2;


u0 = rand(1,1);
r0 = lambda1+((u0.^lambda3-(1-u0).^lambda4)/(-sqrt(lambda2)));
u1 = rand(1,draws);

lambda2t = ones(1,draws).*(omega2+alpha2*r0.^2*lambda2^2+beta2*lambda2);

ret      = ones(1,draws).*(lambda1+(u1.^(lambda3)-(1-u1).^lambda4)/(-sqrt(omega2+alpha2*r0.^2*lambda2t.^2+beta2*lambda2t)));

ewma2    = lambda2t;
ewma3    = zeros(1200,draws);
ewma4    = zeros(1200,draws);
r        = zeros(1001,draws);
e2       = zeros(1001,draws);
e3       = zeros(1001,draws);
e4       = zeros(1001,draws);

for j=1:draws
   for i=1:sampleSize-1
       lambda2t(i+1,j) = omega2+alpha2*ret(i,j)^2*lambda2t(i,j)^2+beta2*lambda2t(i,j);
       u               = rand(1,1);
       ret(i+1,j)      = lambda1+(u.^lambda3-(1-u).^lambda4)/(-sqrt(lambda2t(i+1,j)));
       ewma2(i+1,j)    = (1-theta)*ret(i,j)^2+theta*ewma2(i,j);
       ewma3(i+1,j)    = (1-theta)*ret(i,j)^3+theta*ewma3(i,j);
       ewma4(i+1,j)    = (1-theta)*ret(i,j)^4+theta*ewma4(i,j);
   end
   % burning initial burn observations
   r(:,j)   = ret(burn:sampleSize,j);
   e2(:,j) = ewma2(burn:sampleSize,j);
   e3(:,j) = ewma3(burn:sampleSize,j);
   e4(:,j) = ewma4(burn:sampleSize,j);
end
end



