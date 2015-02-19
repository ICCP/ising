%Plotting routine

load /hfield.out

data=hfield;
x=zeros(1,1000000);
y=zeros(1,1000000);
for i=1:1000000
  x(1,i) = data(i,1);
  y(1,i) = data(i,2);
end

figure(1);
plot(x,y)
title('Magnetization at temp=3')
xlabel('Iterations')
ylabel('Normalized Magnetization')
