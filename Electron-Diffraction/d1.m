%Data
V = [4200 3800 3400 3000 2600 2200 ];     %[V]
R1down = [10.75 11.35 12 12.8 14.15 15.6  ];     %[mm]
R2down = [19.5 20.9 22 23.55 26.25 30 ];   %[mm]
R1upside = [10.5 11.2 12 13 14 15.6 ];     %[mm]
R2upside = [19.5 21 22.5 24.6 27 30];   %[mm]
D = 135; %[mm]
R1=(R1down+R1upside)/2
R2=(R2down+R2upside)/2
d1_theory = 2.13; %[A]
d2_theory = 1.23; %[A]

dr = 0.5; %[mm]
dV = 50; %[V];

%%
%Calculating theory electorn's wave length (de-Brolie)
const = 1.226*10^-9;
lambda_theory = const.*V.^-0.5;
lambda_theory = lambda_theory.*10^10; %[A]

%Calculating measured electorn's wave length (Berg's equation)
lambda1_measured = (d1_theory/D).*R1;
lambda2_measured = (d2_theory/D).*R2;

figure(1)
linfitxy(lambda1_measured, lambda_theory,0.5,0);
title('\lambda theory VS \lambda_1 measured','FontSize',16);
xlabel('\lambda_{measured} [A]','FontSize',14);
ylabel('\lambda_{theory} [A]','FontSize',14);

figure(2)
linfitxy(lambda2_measured, lambda_theory,0.5,0);
title('\lambda theory VS \lambda_2 measured','FontSize',16);
xlabel('\lambda_{measured} [A]','FontSize',14);
ylabel('\lambda_{theory} [A]','FontSize',14);


%%
%Calculating d1 and d2
figure(3)
[p1,g1] = linfitxy(lambda_theory, R1,0,0.5);
title('R_1 VS \lambda theory','FontSize',16);
xlabel('\lambda_{theory} [A]','FontSize',14);
ylabel('R_1 [mm]','FontSize',14);
grid on

figure(4)
[p2,g2] = linfitxy(lambda_theory, R2,0,0.5);
title('R_2 VS \lambda theory','FontSize',16);
xlabel('\lambda_{theory} [A]','FontSize',14);
ylabel('R_2 [mm]','FontSize',14);
grid on

d1_measured = D/p1(1)
d2_measured = D/p2(1)

ed1 = (D^2/p1(1)^4)*g1(1)^2;
ed2 = (D^2/p2(1)^4)*g2(1)^2;

%%
%Calculating Errors
error1 = (1 - d1_measured/d1_theory) * 100
error2 = (1 - d2_measured/d2_theory) * 100

