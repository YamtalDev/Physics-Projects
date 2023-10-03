data1=importdata('run1.txt'); %
data2=importdata('run2.txt'); %
data3=importdata('run3.txt'); %
data4=importdata('run4.txt'); %
v1=data1(:,2); %
v2=data2(:,2); %
v3=data3(:,2); %
v4=data4(:,2); %
I1=data1(:,1); %
I2=data2(:,1); %
I3=data3(:,1); %
I4=data4(:,1); %
figure
plot(v1,I1);
hold on
plot(v2,I2);
hold on
plot(v3,I3);
hold on
plot(v4, I4);
xlabel('Acceleration voltage [Volt] [sec]');%
ylabel('Current  [Amp]');%
title(' Frank Hertz Experiment run:1-4');%


