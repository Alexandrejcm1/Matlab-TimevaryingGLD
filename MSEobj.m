function [SMSE2,SMSE3,SMSE4] = MSEobj(draws,sampleSize,burn,lambda1,lambda4,theta,omega2,alpha2,beta2,lambda3);
draws = 1;
sampleSize=550;
burn  = 50;
lambda1 = 0;
lambda3 = -0.35;
lambda4 = -0.35;
sigma2  = 0.01^2;

% imposing lambda1=0, lambda3=lambda4 => please generalize for any lambda

alpha2 = 0.05;
beta2  = 0.90;
omega2 = 2*(1/(2*lambda3+1)-beta(lambda3+1,lambda4+1))*(1-alpha2-beta2)/sigma2; 
theta = 0.95;


lambda2 = omega2/(1-alpha2-beta2);
u0 = rand(1,draws);
r0 = lambda1+((u0.^lambda3-(1-u0).^lambda4)/(-sqrt(lambda2)));
u1 = rand(1,draws);
% checar contas feitas para chegar ate aqui
lambda2t = ones(1,draws).*(omega2+alpha2*r0.^2+beta2*lambda2);
ret      = ones(1,draws).*(lambda1+(u1.^lambda3-(1-u1).^lambda4)/(-sqrt(omega2+alpha2*r0.^2+beta2*lambda2t)));
ewma2    = lambda2t;
ewma3    = zeros(1,draws);
ewma4    = zeros(1,draws);

for j=1:draws
   for i=1:sampleSize-1
       lambda2t(i+1,j) = omega2+alpha2*ret(i,j)^2+beta2*lambda2t(i,j);
       u               = rand(1,1);
       ret(i+1,j)      = lambda1+(u.^lambda3-(1-u).^lambda4)/(-sqrt(lambda2t(i+1,j)));
       ewma2(i+1,j)    = (1-theta)*ret(i,j)^2+theta*ewma2(i,j);
       ewma3(i+1,j)    = (1-theta)*ret(i,j)^3+theta*ewma3(i,j);
       ewma4(i+1,j)    = (1-theta)*ret(i,j)^4+theta*ewma4(i,j);
end
   % burning initial 50 observations
   r(:,j)   = ret(51:sampleSize,j);
   e2(:,j) = ewma2(51:sampleSize,j);
   e3(:,j) = ewma3(51:sampleSize,j);
   e4(:,j) = ewma4(51:sampleSize,j);
end

% grid for the parameters

%gridlambda3 = [-0.360,-0.355,-0.350,-0.345,-0.340];
%gridomega2 = [0.03,0.04,0.05,0.06,0.07]./100;
%gridalpha2 = [0.045,0.0475,0.05,0.0525,0.055]; 
%gridbeta2  = [0.87,0.885,0.90,0.915,0.93];

% simulate returns and EWMA using grid parameters
J = 100;
%J=1000;
MSE2 = zeros(J,1);
MSE3 = zeros(J,1);
MSE4 = zeros(J,1);


for j=1:1:J
   gridlambda3 = normrnd(-0.35,0.001);
   gridlambda4 = gridlambda3;
   gridomega2 =normrnd(0.0005,0.00001);
   gridalpha2 =normrnd(0.5,0.001);
   gridbeta2  = normrnd(0.90,0.01);

   [rsim,e2sim,e3sim,e4sim] = simulTukey(J,length(r),burn,0,gridlambda4,theta,gridlambda3,gridomega2,gridalpha2,gridbeta2);
   
   % add more unconditional moment conditions
  
    
    MSE2(j,1)          = sum((e2sim(:,:)-e2(:,draws)).^2);
    SMSE2(1,1)         = sum(MSE2(:,1));
    MSE3(1,j)          = sum((e3sim(:,:)-e3(:,draws)).^2);
    SMSE3(1,1)         = sum(MSE3(:,1));
    MSE4(1,j)          = sum((e4sim(:,:)-e4(:,draws)).^2);
    SMSE4(1,1)         = sum(MSE4(:,1));
   %marcelo MSE(j)     = sim((e2sim(:,j)-e2(:,draws))^2
end
end

% matching unconditional moments and MSE's



