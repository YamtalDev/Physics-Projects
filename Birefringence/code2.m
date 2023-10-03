figure
angle=[0,10,20,30,40,50,60,70,80,90];
axisx=(sind(2*angle)).^2;
intensity1=[3.5,3.1,11.5,21.7,23.4,16,8.3,2,0.4,0.3];
plot(axisx,intensity1,'-bo',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.987643 1 .67],...
    'MarkerSize',7)
xlabel('sin^{2}(2\theta)')
ylabel('Measured Light intensity I_f(\theta) [%]')
title('Changing the angle of the poloryser from 2-92 degrees. 4 layers of duct tape. second attempt')
strips=4:10;
errorintensity2=0.2.*[1,1,1,1,1,1,1];
intensity2=[19.6,2.9,25.4, 1.1, 16 ,8, 3.9  ] ; %in percent
figure
plot(strips,intensity2,'k-*')
cftool
xlabel('number of duct tape layers from 4-10')
ylabel('Measured Light intensity I_f(\theta) [%]')
title('Changing the number of layers with a 45 \circ and 135 \circ polarizers ')



