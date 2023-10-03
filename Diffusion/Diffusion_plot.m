Data=load('diff.dat');    %uploading data to matlab
C=Data(:,1)';


c=num2cell(reshape(C, 19*19*19, 2000 ),1); % in c++ we have a grid of 19x19x19, with 2000 steps in t
%so we start sorting by placing each step in time as a seperate cell
% 2000 cells of vectorized 19x19x19 grids, {1:2000}(1x6859) double

for(i=1:length(c))
    SUMPLOT(i)=abs(sum(c{i}));    % vector of the rate of diffusion to compare to analitical
end

%% 2D
for(i=1:length(c))
 C1{i}=reshape(c{i},19, 19, 19);% again reshaping from a 1x6859 vector to 19x19x19 grid
end

%%


for(i=1:length(c)/10) % animation for loop
C1{i*10} = double(squeeze(C1{i*10})); % dimention removal for 'slice' func
h = slice(max(C1{i*10},1.e-15)*100, 1:2:19 , 1:1:19 , 1:2:19); % plot

% basic settings for plot
set(h, 'EdgeColor','none','FaceColor','interp') 
colormap(flipud(jet))
alpha(h,'color') % sets transparacy 
set(gca,'Color','[0.88 0.94 1]')
shading interp %smooth shading
grid off
campos([  (i*100)/length(c) 7 1.4673 ]) %rotation of the video
%campos([ 10  10   180] )            %2d veiw
title('non uniform Rate of diffusion C_{x,y,z,t} with source C_{10,10,10,0}=300  ')
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
 xticklabels({'0', '0.1' , '0.2', '0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
  yticklabels({'0', '0.1' , '0.2', '0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
   zticklabels({'0', '0.25' , '0.5', '0.75','1'})
hold off
pause(0.000000005) % as pause goes to 0 the movie goes faster
end
