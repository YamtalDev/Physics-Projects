
strips=4:10;
errorintensity2=0.2.*[1,1,1,1,1,1,1];
intensity2=[6, 5.6, 4.6, 3.7 ,2.6, 2.1 , 1.4  ] ; %in percent
erroraxis2=0.1.*[1,1,1,1,1,1,1];
figure
plot(strips,intensity2) 
cftool
% [P,SP] =linfitxy(strips,intensity2,erroraxis2,errorintensity2);
xlabel('number of duct tape layers from 4-10')
ylabel('Measured Light intensity I_f(\theta) [%]')
title('Changing the number of layers with a 45 \circ and 135 \circ polarizers ')
figure
angle=[0,10,20,30,40,50,60,70,80,90];
       axisx=(sind(2.*angle)).^2;
       %axisx=cosd(2.*angle);
erroaxisx=0.1.*[1,1,1,1,1,1,1,1,1,1];
errorintensity1=0.2.*[1,1,1,1,1,1,1,1,1,1];
intensity1=[0.1,0.1,0.6,2.3,5.2,9.0,10.3,8.0,4.7,1.0];
plot(axisx,intensity1,'-mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.77 1 .67],...
    'MarkerSize',7)
xlabel('sin^{2}(2\theta)')
ylabel('Measured Light intensity I_f(\theta) [%]')
title('Changing the angle of the poloryser from 2-92 degrees. 4 layers of duct tape.')


