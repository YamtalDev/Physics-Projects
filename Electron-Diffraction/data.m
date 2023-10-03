V = [4200 3800 3400 3000 2600 2200];     %[V]
R1 = [13 13.75 15 15.75 17.5 18.5];     %[mm]
R2 = [23 24.25 25.5 27.75 29.5 32.5];   %[mm]

dr = 0.5; %[mm]
D = 135; %[mm]


%%
const = 1.226*10^-9;
lambda_theory = const.*V.^-0.5;
lambda_theory = lambda_theory.*10^10; %[A]
p1 = linfitxy(lambda_theory, R1,0,0.5);
p2 = linfitxy(lambda_theory, R2,0,0.5);

d_theory = 2.13;
% d1_theory = D/p1(1);
% d2_theory = D/p2(1);



%%
%D = D*10^-3;
%d = d*10^-10;
%R1 = R1.*10^-3;
%R2 = R2.*10^-3;

lambda_measured_R1 = (d_theory/D).*R1;
%lambda_measured_R1 = lambda_measured_R1.*10^10; %[A]

lambda_measured_R2 = (d_theory/D).*R2;
%lambda_measured_R2 = lambda_measured_R2.*10^10; %[A]

%%
dLambda_R1 = abs(lambda_theory - lambda_measured_R1); %[A]
dLambda_R2 = abs(lambda_theory - lambda_measured_R2); %[A]

%%
figure(1)
 linfitxy(lambda_measured_R1,lambda_theory, dLambda_R1, 0);
title('R1');

figure(2)
linfitxy(lambda_measured_R2,lambda_theory, dLambda_R2, 0);
title('R2');

%%
g1 = linfitxy(R1,lambda_theory,0.5,0);
g2 = linfitxy(R2,lambda_theory,0.5,0);

%%
d1_measured = D/g1(1);
d2_measured = D/g2(1);