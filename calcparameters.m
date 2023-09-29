clear all;
clc;


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

J=100;
MSE2 = zeros(J,1);
MSE3 = zeros(J,1);
MSE4 = zeros(J,1);


% matching unconditional moments and MSE's
%x0 = [1000,0.05,0.9,-0.35];
%[omega2hat,alpha2hat,beta2hat,lambda3hat] = fminsearch(@MSEobj,arguments,options)
% tudo que foi necessario p rodar a mse menos os parametros q quero
%[omega2hat,alpha2hat,beta2hat,lambda3hat] = fminsearch(@MSEobj,x0,J,length(r),burn,0,gridlambda4,theta,options);



options=optimset('MaxIter', 100, 'MaxFunEvals', 100);
x0 = [0.0005,0.05,0.9,-0.35];
%[omega2hat,alpha2hat,beta2hat,lambda3hat] = fminsearch(@MSEobj,x0,draws,sampleSize,burn,lambda1,lambda4,theta,options);

%[omega2hat,alpha2hat,beta2hat,lambda3hat] = fminsearch(@MSEobj,arguments,options)
% acho que isso ta certo, agora fazer rodar isso
%[omega2hat,alpha2hat,beta2hat,lambda3hat] = fminsearch(@MSEobj,x0,draws,sampleSize,burn,lambda1,lambda4,theta,options);
%estimou, o problema ta na matriz J, testando p3 MSE
[x,fval,exitflag,output] = fminsearch(@MSEobj,x0,draws,sampleSize,burn,lambda1,lambda4,theta,options)


omega2hat=x(1)
alpha2hat=x(2)
beta2hat=x(3)
lambda3hat=x(4)

