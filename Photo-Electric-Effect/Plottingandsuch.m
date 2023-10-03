Datap=load('Navier_stokes.dat');    
Datau=load('Navier_stokes2.dat');    
Datav=load('Navier_stokes3.dat');   
%%
n=55
k=length(Datap)/(n*n)
cp=num2cell(reshape(Datap, n*n, k ),1);
cu=num2cell(reshape(Datau, n*n, k ),1);
cv=num2cell(reshape(Datav, n*n, k ),1);
for(i=1:length(cp))
 C1p{i}=reshape(cp{i},n, n);% again reshaping from a 1x6859 vector to 19x19x19 grid
end
for(i=1:length(cp))
 C1u{i}=reshape(cu{i},n, n);% again reshaping from a 1x6859 vector to 19x19x19 grid
end
for(i=1:length(cp))
 C1v{i}=reshape(cv{i},n, n);% again reshaping from a 1x6859 vector to 19x19x19 grid
end


%%

for i=1:length(C1p)
   %contourf(Datap([n*i:n*(i+1)],[1:n]))
    hold on
    contourf(C1p{i})
    quiver(C1u{i},C1v{i},'k')
    set(gcf, 'WindowState', 'maximized');
    hold off
pause(0.00005) 

end