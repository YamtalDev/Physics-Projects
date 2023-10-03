Datap=load('Navier_stokes.dat');    %Loading In Data that was processed by Matlab
Datau=load('Navier_stokes2.dat');    % "
Datav=load('Navier_stokes3.dat');    % "
%%
n=199;                               % Gridsize from c++ - 1
G=length(Datav)/(n*n);               % Finding how many instances we have
%%
cp=num2cell(reshape(Datap, n*n, G ),1);  %rescaling/orgainzing vectors for plotting

cu=num2cell(reshape(Datau, n*n, G ),1); %rescaling/orgainzing vectors for plotting
 
cv=num2cell(reshape(Datav, n*n, G ),1); %rescaling/orgainzing vectors for plotting
%%
for(i=1:length(cp))
 C1p{i}=reshape(cp{i},n, n);% again reshaping from a 1xn^2 vector to nxn grid
end
for(i=1:length(cp))
 C1u{i}=reshape(cu{i},n, n);% again reshaping from a 1xn^2 vector to nxn grid
end
for(i=1:length(cp))
 C1v{i}=reshape(cv{i},n, n);% again reshaping from a 1xn^2 vector to nxn grid
end
%% VIDEO PRINTING
% the videos where quite large, so we must split the videos in four parts, 
% j for 1:G, i for video purposes. 

j=1;                              
for i=1:1500
  hold on

     %contourf(C1p{j},x,':'); 
     q=pcolor(C1u{j}+C1v{j});
     set(q, 'EdgeColor', 'none');
     %h=quiver(C1u{j},C1v{j},5);
     %set( h, 'Color', 'w' )
    axis([0 200 0 200])
     xlabel('i')
     ylabel('j')
     
     c = colorbar;
     c.Label.String = ' U_x Velocity';
     %caxis([0 0.7]) 
     %set(gcf, 'Position', get(0, 'Screensize'));
     drawnow
     frame(i)=getframe(gcf);
     hold off
     cla reset; 
     j=j+1
end
%%
video = VideoWriter( 'part1.avi' , 'Motion JPEG AVI' ) ; 
video.FrameRate=30;
open(video)
writeVideo(video, frame);
close(video)
%%
j=1500;
for i=1:1500
  hold on
  
     %contourf(C1p{j},x,':'); 
     q=pcolor(C1u{j}+C1v{j});
     set(q, 'EdgeColor', 'none');
     %h=quiver(C1u{j},C1v{j},4);
     %set( h, 'Color', 'w' )
    axis([0 200 0 200])
     xlabel('i')
     ylabel('j')
     
     c = colorbar;
     c.Label.String = ' U_x Velocity';
     %caxis([min(Datap) max(Datap)]) 
     %set(gcf, 'Position', get(0, 'Screensize'));
     drawnow
     frame(i)=getframe(gcf);
     hold off
     cla reset; 
     j=j+1
end
%%
video = VideoWriter( 'part2.avi' , 'Motion JPEG AVI' ) ; 
video.FrameRate=30;
open(video)
writeVideo(video, frame);
close(video)
pause(15)
%%
j=3000;
for i=1:1500
  hold on
  
     %contourf(C1p{j},x,':'); 
     q=pcolor(C1u{j}+C1v{j});
     set(q, 'EdgeColor', 'none');
     h=quiver(C1u{j},C1v{j},5);
     set( h, 'Color', 'w' )
    axis([0 200 0 200])
     xlabel('i')
     ylabel('j')
     
     c = colorbar;
     c.Label.String = ' U_x Velocity';
     %caxis([min(Datap) max(Datap)]) 
     %set(gcf, 'Position', get(0, 'Screensize'));
     drawnow
     frame(i)=getframe(gcf);
     hold off
     cla reset; 
     j=j+1
end
%%
video = VideoWriter( 'part3.avi' , 'Motion JPEG AVI' ) ; 
video.FrameRate=30;
open(video)
writeVideo(video, frame);
close(video)
pause(15)
%%
j=4500;
for i=1:1500
  hold on
  
     %contourf(C1p{j},x,':'); 
     q=pcolor(C1u{j}+C1v{j});
     set(q, 'EdgeColor', 'none');
     h=quiver(C1u{j},C1v{j},7);
     set( h, 'Color', 'w' )
    axis([0 200 0 200])
     xlabel('i')
     ylabel('j')
     
     c = colorbar;
     c.Label.String = ' U_x Velocity';
     %caxis([min(Datap) max(Datap)]) 
     %set(gcf, 'Position', get(0, 'Screensize'));
     drawnow
     frame(i)=getframe(gcf);
     hold off
     cla reset; 
     j=j+1
end
video = VideoWriter( 'part4.avi' , 'Motion JPEG AVI' ) ; 
video.FrameRate=30;
open(video)
writeVideo(video, frame);
close(video)
pause(15)
%% Plots
U=C1u{end};
V=C1v{end};
P=C1p{end};
figure(7)
contourf(P);
hold on 
 o=streamslice(U,V,2);

  set( o, 'Color', 'w' )
       axis([0 150 0 150])
     xlabel('i')
     ylabel('j')
   
     
     c = colorbar;
     c.Label.String = 'Pressure';
title('Final State T_{final} Re=100k')
figure(2)
pcolor(U+V);
hold on 
 o=streamslice(U,V,2);

  set( o, 'Color', 'w' )
       axis([0 150 0 150])
     xlabel('i')
     ylabel('j')
   
     
     c = colorbar;
     c.Label.String = 'magnitube of velocity';
title('Final State T_{final}  Re=100k')

figure(3)
centerline=U(:,75)+V(:,75)
 plot(centerline,[1:199])
  title('U (flow speed) at i=i/2 Re=100k')
  xlabel('U(Flow Speed)')   
  ylabel('j')
  
