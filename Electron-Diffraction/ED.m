b=8.2e-03;
d=8e-03;
dd=0.5e-03
g=9.80665;
p=10^5;
V=398.5;
dv=0.2;
n=18.42e-06;
dn=0.001e-05;
rho=886;
q0=1.602e-19;
dis=3*10^(-4);
E=V/d;

tr = [0.93,0.95,0.94,0.9,0.92];
tf = [2.63,2.62,2.64,2.6,2.66];

tf2 = [2.52,2.48,2.55,2.48,2.59];
tr2 = [2.62,2.43,2.53,2.5,2.57];

tf3 = [2.61,2.5,2.51,2.45,2.5];
tr3 = [2,1.9,2.1,2.07,2.11];

tf4 = [2.4,2.4,2.4,2.22,2.3];
tr4 = [1.58,1.49,1.57,1.55,1.52];

tf5 = [1.45,1.53,1.65,1.61,1.5];
tr5 = [2.65,2.54,2.47,2.336,2.44];


tfm = sum(tf)/length(tf);
trm = sum(tr)/length(tr);


tfm2 = sum(tf2)/length(tf2);
trm2 = sum(tr2)/length(tr2);


tfm3 = sum(tf3)/length(tf3);
trm3 = sum(tr3)/length(tr3);


tfm4 = sum(tf4)/length(tf4);
trm4 = sum(tr4)/length(tr4);


tfm5 = sum(tf5)/length(tf5);
trm5 = sum(tr5)/length(tr5);


%% 
Vf = dis/tfm;
Vr = dis/trm;

Vf2 = dis/tfm2;
Vr2 = dis/trm2;

Vf3 = dis/tfm3;
Vr3 = dis/trm3;

Vf4 = dis/tfm4;
Vr4 = dis/trm4;

Vf5 = dis/tfm5;
Vr5 = dis/trm5;

a=sqrt((9*n*Vf)/(2*rho*g));
a2=sqrt((9*n*Vf2)/(2*rho*g));
a3=sqrt((9*n*Vf3)/(2*rho*g));
a4=sqrt((9*n*Vf4)/(2*rho*g));
a5=sqrt((9*n*Vf5)/(2*rho*g));

ac=a*sqrt((p*a)/(p*a+b));
ac2=a2*sqrt((p*a2)/(p*a2+b));
ac3=a3*sqrt((p*a3)/(p*a3+b));
ac4=a4*sqrt((p*a4)/(p*a4+b));
ac5=a5*sqrt((p*a5)/(p*a5+b));

q_radius=(ac^3)*((4*pi*rho*g)/3)*((Vf+Vr)/(E*Vf));
q_radius2=(ac2^3)*((4*pi*rho*g)/3)*((Vf2+Vr2)/(E*Vf2));
q_radius3=(ac3^3)*((4*pi*rho*g)/3)*((Vf3+Vr3)/(E*Vf3));
q_radius4=(ac4^3)*((4*pi*rho*g)/3)*((Vf4+Vr4)/(E*Vf4));
q_radius5=(ac5^3)*((4*pi*rho*g)/3)*((Vf5+Vr5)/(E*Vf5));

Qnum_radius=[q_radius/q0,q_radius2/q0,q_radius3/q0,q_radius4/q0,q_radius5/q0];


error = 100.*(abs(1 - Qnum_radius))

q=((4*pi/3)*(sqrt((b/(2*p))^2+(9*n*Vf)/(2*rho*g))-b/(2*p))^3)*(((rho*g*d)*(Vf+Vr))/(V*Vf));
q2=((4*pi/3)*(sqrt((b/(2*p))^2+(9*n*Vf2)/(2*rho*g))-b/(2*p))^3)*(((rho*g*d)*(Vf2+Vr2))/(V*Vf2));
q3=((4*pi/3)*(sqrt((b/(2*p))^2+(9*n*Vf3)/(2*rho*g))-b/(2*p))^3)*(((rho*g*d)*(Vf3+Vr3))/(V*Vf3));
q4=((4*pi/3)*(sqrt((b/(2*p))^2+(9*n*Vf4)/(2*rho*g))-b/(2*p))^3)*(((rho*g*d)*(Vf4+Vr4))/(V*Vf4));
q5=((4*pi/3)*(sqrt((b/(2*p))^2+(9*n*Vf5)/(2*rho*g))-b/(2*p))^3)*(((rho*g*d)*(Vf5+Vr5))/(V*Vf5));

Qnum=[q/q0,q2/q0,q3/q0,q4/q0,q5/q0]
error = 100.*(abs(1 - Qnum))