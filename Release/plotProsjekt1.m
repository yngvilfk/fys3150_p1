close all
test10

figure(1)
subplot(3,1,1)
plot(x_vec,u_vec,'-r', 'LineWidth',2)
hold on
plot(x_vec,v_vec,'ob', 'LineWidth',2)
xlabel('x')
ylabel('u(x)')
legend('analytical solution','numerical solution')
title(sprintf('n = %s',num2str(length(x_vec)-2)))
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 100)


test100


subplot(3,1,2)
plot(x_vec,u_vec,'-r', 'LineWidth',2)
hold on
plot(x_vec,v_vec,'ob', 'LineWidth',2)
xlabel('x')
ylabel('u(x)')
legend('analytical solution','numerical solution')
title(sprintf('n = %s',num2str(length(x_vec)-2)))
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 100)


test1000


subplot(3,1,3)
plot(x_vec,u_vec,'-r', 'LineWidth',2)
hold on
plot(x_vec,v_vec,'ob', 'LineWidth',2)
xlabel('x')
ylabel('u(x)')
legend('analytical solution','numerical solution')
title(sprintf('n = %s',num2str(length(x_vec)-2)))
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 100)

print -dpng  'project1.png' -r100

