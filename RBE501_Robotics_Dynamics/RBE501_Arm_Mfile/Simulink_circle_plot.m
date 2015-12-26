figure(5) 
plot3(Pee_circle.Data(:,1),Pee_circle.Data(:,3),Pee_circle.Data(:,2),'LineWidth',3)
%axis([0,0.2,0,0.3,0.2,0.3])
title('Simulink End-effector Position (circle)','FontSize',15)
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on