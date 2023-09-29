tic
clear all

% initial parameters values

theta   = 0.95;
alpha2  = 0.05;
beta2   = 0.90;
sigma2  = 0.01^2;
lambda3 = -0.15;
lambda4 = -0.2;

A = 1/(lambda3+1)-1/(lambda4+1);
B = 1/(2*lambda3+1)+1/(2*lambda4+1)-2*beta(lambda3+1,lambda4+1);
C = 1/(3*lambda3+1)-1/(3*lambda4+1)-3*beta(2*lambda3+1,lambda4+1)+3*beta(lambda3+1,2*lambda4+1);
D = 1/(4*lambda3+1)+1/(4*lambda4+1)-4*beta(3*lambda3+1,lambda4+1)+6*beta(2*lambda3+1,2*lambda4+1)-4*beta(lambda3+1,3*lambda4+1);

lambda2 = sqrt((B-A^2)/sigma2);
omega2  = (1-beta2)*lambda2-alpha2*(B-A^2);
lambda1 = -A/lambda2;

% simulating returns

draws      = 2000;   % number of MC replications, at least 2000!!!
sampleSize = 1200;  % try 1200, burn the first (700 or 200) to end up with sample sizes of (500,1000)
burn       = 700;

%u0 = rand(1,draws);
u0 = rand(1,1);
r0 = lambda1+((u0.^lambda3-(1-u0).^lambda4)/(-sqrt(lambda2)));
u1 = rand(1,draws);
%Aqui nao eh raiz - lambda2 mas sim lambda2 pq nao eh variancia mas sim
%lambda2 somente ( testar)

lambda2t = ones(1,draws).*(omega2+alpha2*r0.^2*lambda2.^2+beta2*lambda2);
%ret      = ones(1,draws).*(lambda1+(u1.^lambda3-(1-u1).^lambda4)/(-sqrt(omega2+alpha2*r0.^2*lambda2t.^2+beta2*lambda2t)));
ret      = ones(1,draws).*(lambda1+(u1.^lambda3-(1-u1).^lambda4)/(omega2+alpha2*r0.^2*lambda2t.^2+beta2*lambda2t));

ewma2    = zeros(1,draws);
ewma3    = zeros(1,draws);
ewma4    = zeros(1,draws);

for j=1:draws
   for i=1:sampleSize-1
       lambda2t(i+1,j) = omega2+alpha2*ret(i,j)^2*lambda2t(i,j)^2+beta2*lambda2t(i,j);
       u               = rand(1,1);
  %    ret(i+1,j)      = lambda1+(u.^lambda3-(1-u).^lambda4)/(-sqrt(lambda2t(i+1,j)));
       ret(i+1,j)      = lambda1+(u.^lambda3-(1-u).^lambda4)/(lambda2t(i+1,j));

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

% indirect inference
% simulate returns and EWMA using grid parameters

J    = 1000;         % grid size (try also 1000)
MSE2 = zeros(J,1);
MSE3 = zeros(J,1);
MSE4 = zeros(J,1);

for j=1:1:J
   gridlambda3 = -betarnd(3,3)*0.25;
   histlambda3(j,:) = gridlambda3;
   gridlambda4 = -betarnd(3,3)*0.25;
   histlambda4(j,:) = gridlambda4;
   gridsigma2  = normrnd(0.001,0.001);
   histsigma2(j,:) = gridsigma2;
   gridalpha2  = normrnd(0.05,0.01);
   histalpha2(j,:) = gridalpha2;
   gridbeta2 = betarnd(25,1.5)-gridalpha2; 
   histbeta2(j,:) = gridbeta2;

   [rsim,e2sim,e3sim,e4sim] = simulTukey(sampleSize,burn,gridlambda4,theta,gridlambda3,gridsigma2,gridalpha2,gridbeta2);
   
   % MSE's   (antes pegava somente a primeira coluna , agora a matriz toda( vai demorar mais))
   
   MSE2(j,:) = immse(e2sim(:,:),e2(:,:));
   MSE3(j,:) = immse(e3sim(:,:),e3(:,:));
   MSE4(j,:) = immse(e4sim(:,:),e4(:,:));
   %MSE2(j,1) = mean((e2sim(:,1)-e2(:,draws)).^2);
   %MSE3(j,1) = mean((e3sim(:,1)-e3(:,draws)).^2);
   %MSE4(j,1) = mean((e4sim(:,1)-e4(:,draws)).^2);
   MMSE(j,1) = (MSE2(j,1)+MSE3(j,1)+MSE4(j,1))/3;
   end

[Mhat,Ihat] = min(MSE2);
[Mtilde,Itilde] = min(MMSE);

alpha2hat  = histalpha2(Ihat,:);
beta2hat   = histbeta2(Ihat,:);
lambda3hat = histlambda3(Ihat,:);
sigma2hat  = histsigma2(Ihat,:);
lambda4hat = histlambda4(Ihat,:);

Ahat    = 1/(lambda3hat+1)-1/(lambda4hat+1);
Bhat    = 1/(2*lambda3hat+1)+1/(2*lambda4hat+1)-2*beta(lambda3hat+1,lambda4hat+1);
Chat    = 1/(3*lambda3hat+1)-1/(3*lambda4hat+1)-3*beta(2*lambda3hat+1,lambda4hat+1)+3*beta(lambda3hat+1,2*lambda4hat+1);
Dhat    = 1/(4*lambda3hat+1)+1/(4*lambda4hat+1)-4*beta(3*lambda3hat+1,lambda4hat+1)+6*beta(2*lambda3hat+1,2*lambda4hat+1)-4*beta(lambda3hat+1,3*lambda4hat+1);
skewhat = (1/(Bhat-Ahat^2)^1.5)*(Chat-3*Ahat*Bhat-2*Ahat^3);
kurthat = ((1/Bhat-Ahat^2)^2)*(Dhat-4*Ahat*Chat+6*Ahat^2*Bhat+3*Ahat^4);

yhat = [alpha2hat,beta2hat,sigma2hat,lambda3hat,lambda4hat,skewhat,kurthat,Mhat]

alpha2tilde  = histalpha2(Itilde,:);
beta2tilde   = histbeta2(Itilde,:);
lambda3tilde = histlambda3(Itilde,:);
sigma2tilde  = histsigma2(Itilde,:);
lambda4tilde = histlambda4(Itilde,:);

Atilde    = 1/(lambda3tilde+1)-1/(lambda4tilde+1);
Btilde    = 1/(2*lambda3tilde+1)+1/(2*lambda4tilde+1)-2*beta(lambda3tilde+1,lambda4tilde+1);
Ctilde    = 1/(3*lambda3tilde+1)-1/(3*lambda4tilde+1)-3*beta(2*lambda3tilde+1,lambda4tilde+1)+3*beta(lambda3tilde+1,2*lambda4tilde+1);
Dtilde    = 1/(4*lambda3tilde+1)+1/(4*lambda4tilde+1)-4*beta(3*lambda3tilde+1,lambda4tilde+1)+6*beta(2*lambda3tilde+1,2*lambda4tilde+1)-4*beta(lambda3tilde+1,3*lambda4tilde+1);
skewtilde = (1/(Btilde-Atilde^2)^1.5)*(Ctilde-3*Atilde*Btilde-2*Atilde^3);
kurttilde = ((1/Btilde-Atilde^2)^2)*(Dtilde-4*Atilde*Ctilde+6*Atilde^2*Btilde+3*Atilde^4);

ytilde  = [alpha2tilde,beta2tilde,sigma2tilde,lambda3tilde,lambda4tilde,skewtilde,kurttilde,Mtilde]

% incluir calculo de vies e erro quadratico medio para cada parametro

alpha2hatbias    = alpha2hat   - alpha2;
beta2hatbias     = beta2hat    - beta2;
lambda3hatbias   = lambda3hat  - lambda3;
sigma2hatbias    = sigma2hat   - sigma2;
lambda4hatbias   = lambda4hat  - lambda4;

yhatbias = [alpha2hatbias,beta2hatbias,lambda3hatbias,sigma2hatbias,lambda4hatbias]

alpha2tildebias  = alpha2tilde   - alpha2;
beta2tildebias   = beta2tilde    - beta2;
lambda3tildebias = lambda3tilde  - lambda3;
sigma2tildebias  = sigma2tilde   - sigma2;
lambda4tildebias = lambda4tilde - lambda4;

ytildebias = [alpha2tildebias,beta2tildebias,lambda3tildebias,sigma2tildebias,lambda4tildebias]

MSEalpha2hat     =  var(histalpha2)  -  alpha2hatbias^2; 
MSEbeta2hat      =  var(histbeta2)   -  beta2hatbias^2; 
MSElambda3hat    =  var(histlambda3) -  lambda3hatbias^2;
MSEsigma2hat     =  var(histsigma2)  -  sigma2hatbias^2;
MSElambda4hat    =  var(histlambda4) -  lambda4hatbias^2;

MSEyhat = [MSEalpha2hat,MSEbeta2hat,MSElambda3hat,MSEsigma2hat,MSElambda4hat]

MSEalpha2tilde    =  var(histalpha2)  -  alpha2tildebias^2;
MSEbeta2tilde     =  var(histbeta2)   -  beta2tildebias^2;
MSElambda3tilde   =  var(histlambda3) -  lambda3tildebias^2;
MSEsigma2tilde    =  var(histsigma2)  -  sigma2tildebias^2;
MSElambda4tilde   =  var(histlambda4) -  lambda4tildebias^2;

MSEytilde = [MSEalpha2tilde,MSEbeta2tilde,MSElambda3tilde,MSEsigma2tilde,MSElambda4tilde]

resulthat = [alpha2hat,alpha2hatbias,MSEalpha2hat, beta2hat,beta2hatbias,MSEbeta2hat,sigma2hat,sigma2hatbias,MSEsigma2hat,lambda3hat,lambda3hatbias,MSElambda3hat,lambda4hat,lambda4hatbias,MSElambda4hat,skewhat,kurthat,Mhat]
resulttilde = [alpha2tilde,alpha2tildebias,MSEalpha2tilde, beta2tilde,beta2tildebias,MSEbeta2tilde,sigma2tilde,sigma2tildebias,MSEsigma2tilde,lambda3tilde,lambda3tildebias,MSElambda3tilde,lambda4tilde,lambda4tildebias,MSElambda4tilde,skewtilde,kurttilde,Mtilde]

toc