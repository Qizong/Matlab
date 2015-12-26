%Model the Serpenoid curve that
%---- Use Serpentine locomotion
%---- Depend on a, b and c parameters
%---- by Qizong Wu
clc
clear all
% Snake Locomotion Parameters
n =20;         % number of segments
L = 32;         % lenght of the snake (not used for now)
a1 = pi/3;      % degree of undulation
a2 = pi/2;
a3 = pi*2/3;
a1DEG = a1*(180/pi);
a2DEG = a2*(180/pi);
a3DEG = a3*(180/pi);
b1 = 2*pi;       % periods of unit lenghts
b2 = 4*pi;
b3 = 6*pi;
c1 = 0;          % changes the direction of propagation
c2 = pi/2;
c3 = pi;
omega=.5        % in rad/s
frequency= omega/(2*pi) %in Hz
lag = (2*pi)/n	% in rad
% The X-Y position of the center of each segment can
% be aproximated by the following code:
x1(1) = 0;      %For the origin or tail position
y1(1) = 0;      %For the origin or tail position
x2(1) = 0;      
y2(1) = 0;      
x3(1) = 0;      
y3(1) = 0;      
%position is plotted in the graph
%to calculate the remaining x-y positions following the formulas
figure(1)
% a varies
subplot(3,1,1)
for k=1:n
    x1(k+1) = x1(k)+(1/n)*cos(a1*cos(k*b1/n)); %x position for i-th segment
    y1(k+1) = y1(k)+(1/n)*sin(a1*cos(k*b1/n)); %y position for i-th segment
    x2(k+1) = x2(k)+(1/n)*cos(a2*cos(k*b1/n)); %x position for i-th segment
    y2(k+1) = y2(k)+(1/n)*sin(a2*cos(k*b1/n)); %y position for i-th segment
    x3(k+1) = x3(k)+(1/n)*cos(a3*cos(k*b1/n)); %x position for i-th segment
    y3(k+1) = y3(k)+(1/n)*sin(a3*cos(k*b1/n)); %y position for i-th segment
end
plot(x1,y1,'b-',x2,y2,'b--',x3,y3,'b:'); %plot the x-y coordinates
xlabel('x')
ylabel('y')
title('Postion Graph in Cartesian Plane (b = 2\pi; c = 0)')
legend('a = \pi/3','a = \pi/2','a = 2\pi/3')
% b varies
subplot(3,1,2)
for k=1:n
    x1(k+1) = x1(k)+(1/n)*cos(a1*cos(k*b1/n)); %x position for i-th segment
    y1(k+1) = y1(k)+(1/n)*sin(a1*cos(k*b1/n)); %y position for i-th segment
    x2(k+1) = x2(k)+(1/n)*cos(a1*cos(k*b2/n)); %x position for i-th segment
    y2(k+1) = y2(k)+(1/n)*sin(a1*cos(k*b2/n)); %y position for i-th segment
    x3(k+1) = x3(k)+(1/n)*cos(a1*cos(k*b3/n)); %x position for i-th segment
    y3(k+1) = y3(k)+(1/n)*sin(a1*cos(k*b3/n)); %y position for i-th segment
end
plot(x1,y1,'b-',x2,y2,'b--',x3,y3,'b:'); %plot the x-y coordinates
xlabel('x')
ylabel('y')
title('Postion Graph in Cartesian Plane (a = \pi/3; c = 0)')
legend('b = 2\pi','b = 4\pi','b = 6\pi')
% c varies
subplot(3,1,3)
for k=1:n
    x1(k+1) = x1(k)+(1/n)*cos(a1*cos(k*b3/n)+k*c1/n); %x position for i-th segment
    y1(k+1) = y1(k)+(1/n)*sin(a1*cos(k*b3/n)+k*c1/n); %y position for i-th segment
    x2(k+1) = x2(k)+(1/n)*cos(a1*cos(k*b3/n)+k*c2/n); %x position for i-th segment
    y2(k+1) = y2(k)+(1/n)*sin(a1*cos(k*b3/n)+k*c2/n); %y position for i-th segment
    x3(k+1) = x3(k)+(1/n)*cos(a1*cos(k*b3/n)+k*c3/n); %x position for i-th segment
    y3(k+1) = y3(k)+(1/n)*sin(a1*cos(k*b3/n)+k*c3/n); %y position for i-th segment
end
plot(x1,y1,'b-',x2,y2,'b--',x3,y3,'b:'); %plot the x-y coordinates
xlabel('x')
ylabel('y')
title('Postion Graph in Cartesian Plane (a = \pi/3; b = 6\pi)')
legend('c = 0','c = pi/2','c = \pi')

%% 8 Segments Snake Model
clc
clear
hold on
n =8                %number of segments
L = 32;             %lenght of the snake (not used for now)
a = pi/3            %degree of undulation
aDEG = a*(180/pi)
b = 2*pi            %periods of unit lenghts
c = 0;              %changes the direction of propagation
omega=.5            %in rad/s
frequency= omega/(2*pi) %in Hz
lag = (2*pi)/n      %in rad
% The X-Y position of the center of each segment can
% be aproximated by the following code:
x(1) = 0;           %For the origin or tail position
y(1) = 0;           %For the origin or tail position
%position is plotted in the graph
%to calculate the remaining x-y positions following the formulas
%in the paper "Snake Analysis"
figure(2)
for i=1:1:n+1
    for k=1:1:i
        x(k+1) = x(k)+(1/n)*cos(a*cos(k*b/n)); %x position for i-th segment
        y(k+1) = y(k)+(1/n)*sin(a*cos(k*b/n)); %y position for i-th segment
        if i<=n
        line([x(k),x(k+1)],[y(k),y(k+1)],'LineWidth',3)
        end
    end
hold on
end
plot(x(2:8),y(2:8),'ro','MarkerFaceColor','r'); %plot the x-y coordinates
title('Snake Model with 8 Segments (a = \pi/3; b = 2\pi; c = 0)')
%pause(.05); %pause to see it like an animation
hold off
%% Calculate the change in theta with time
% Theta angles for each joint changing with time for 20 seconds
figure(3)
for t=1:1:20
for i=1:1:n
theta(t,i) = 2*a*(180/pi)*abs(sin(b/(2*n)))*sin(omega*t+(b/n)*(i-1));
end
%Plot the theta angle with time
hold on
plot(t,theta(t,1),'-or')
plot(t,theta(t,2),'-ob')
plot(t,theta(t,3),'-om')
plot(t,theta(t,4),'-oy')
plot(t,theta(t,5),'-og')
plot(t,theta(t,6),'-ok')
plot(t,theta(t,7),'-oc')
xlabel('Time (s)')
ylabel('\theta [degrees]')
title('Change of Angle \theta (n = 8)')
end
for i=1:1:n-1
line([0,1],[0,theta(1,i)],'color','k')
for t=1:1:19
line([t,t+1],[theta(t,i), theta(t+1,i)],'color','k')
end
end
hold off

%% 6 plots of Snake motion 
x(1) = 0; %For the origin or tail position
y(1) = 0; %For the origin or tail position
%position is plotted in the graph
figure(4)
%to calculate the remaining x-y positions following the formulas
%in the paper "Snake Analysis"
for t=1:2:11
%figure
subplot(3,2,(t+1)/2)
    for i=1:1:n
        for k=1:1:i
            x(k+1) = x(k)+(1/n)*cos(theta(t,k)*(pi/180)); %x position for i-th segment
            y(k+1) = y(k)+(1/n)*sin(theta(t,k)*(pi/180)); %y position for i-th segment
                if i<=n
                    line([x(k),x(k+1)],[y(k),y(k+1)],'LineWidth',2)
                end
        end
        hold on
        plot(x(k),y(k),'ro','MarkerFaceColor','r'); %plot the x-y coordinates
        xlabel('x ')
        ylabel('y ')
        title(['Postion (t = ',num2str(t),' s)'])
        axis([0 1 -.3 .3])
        grid on
    end 
    mov(t) = getframe;
end
%% Snake motion Movie
x(1) = 0; %For the origin or tail position
y(1) = 0; %For the origin or tail position
%plot(x(1),y(1),'bo') %plot the tail position in the graph
hold on %hold this position in the graph when another point
%position is plotted in the graph

%to calculate the remaining x-y positions following the formulas
%in the paper "Snake Analysis"
for t=1:20
figure

    for i=1:1:n
        for k=1:1:i
            x(k+1) = x(k)+(1/n)*cos(theta(t,k)*(pi/180)); %x position for i-th segment
            y(k+1) = y(k)+(1/n)*sin(theta(t,k)*(pi/180)); %y position for i-th segment
                if i<=n
                    line([x(k),x(k+1)],[y(k),y(k+1)],'LineWidth',2)
                end
        end
        hold on
        plot(x(k),y(k),'ro','MarkerFaceColor','r'); %plot the x-y coordinates
        xlabel('x ')
        ylabel('y ')
        title('Shape of a 8-link Snake Robot')
        axis([0 1 -.3 .3])
        grid on
    end 
    mov(t) = getframe;
end
% Moving animation
figure(33)
pause(2);
movie(mov,1,1)
% %to calculate the delay in seconds of the control signals
% % 'Time delay between square-wave control signals..........................'
delay = 2*pi*(b/(n*2*pi))/omega
%% Angular Velocity Calculation
D = [1 -1 0 0 0 0 0 0;
     0 1 -1 0 0 0 0 0;
     0 0 1 -1 0 0 0 0;
     0 0 0 1 -1 0 0 0;
     0 0 0 0 1 -1 0 0;
     0 0 0 0 0 1 -1 0;
     0 0 0 0 0 0 1 -1];
L = 0.32;
tt = 20;
for i = 1:tt
    L_m(:,:,i) = [L/2 0 0 0 0 0 0 0;
              L*cos(theta(i,2)-theta(i,1)) L/2 0 0 0 0 0 0;
              L*cos(theta(i,3)-theta(i,1)) L*cos(theta(i,3)-theta(i,2)) L/2 0 0 0 0 0;
              L*cos(theta(i,4)-theta(i,1)) L*cos(theta(i,4)-theta(i,2)) L*cos(theta(i,4)-theta(i,3)) L/2 0 0 0 0;
              L*cos(theta(i,5)-theta(i,1)) L*cos(theta(i,5)-theta(i,2)) L*cos(theta(i,5)-theta(i,3)) L*cos(theta(i,5)-theta(i,4)) L/2 0 0 0;
              L*cos(theta(i,6)-theta(i,1)) L*cos(theta(i,6)-theta(i,2)) L*cos(theta(i,6)-theta(i,3)) L*cos(theta(i,6)-theta(i,4)) L*cos(theta(i,6)-theta(i,5)) L/2 0 0;
              L*cos(theta(i,7)-theta(i,1)) L*cos(theta(i,7)-theta(i,2)) L*cos(theta(i,7)-theta(i,3)) L*cos(theta(i,7)-theta(i,4)) L*cos(theta(i,7)-theta(i,5)) L*cos(theta(i,7)-theta(i,6)) L/2 0;
              L*cos(theta(i,8)-theta(i,1)) L*cos(theta(i,8)-theta(i,2)) L*cos(theta(i,8)-theta(i,3)) L*cos(theta(i,8)-theta(i,4)) L*cos(theta(i,8)-theta(i,5)) L*cos(theta(i,8)-theta(i,6)) L*cos(theta(i,8)-theta(i,7)) L/2]
    SC(:,:,i) = [-sin(theta(i,1)) cos(theta(i,1));
          -sin(theta(i,2)) cos(theta(i,2));
          -sin(theta(i,3)) cos(theta(i,3));
          -sin(theta(i,4)) cos(theta(i,4));
          -sin(theta(i,5)) cos(theta(i,5));
          -sin(theta(i,6)) cos(theta(i,6));
          -sin(theta(i,7)) cos(theta(i,7));
          -sin(theta(i,8)) cos(theta(i,8))];
    F(:,:,i) = D * L_m(:,:,i) * SC(:,:,i);
    phid(:,i) = F(:,:,i) * [2;0];
end
figure(6)
subplot(3,3,1)
xt = 1:20;
yt = phid(1,:);
xt1 = 1:0.1:20;
yt1 = spline(xt,yt,xt1);
plot(xt1,yt1)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
title({'Joint 1 Angular Velocity';'(desired speed Vx = 2m/s, Vy = 0)'})

subplot(3,3,2)
xt = 1:20;
yt = phid(2,:);
xt1 = 1:0.1:20;
yt1 = spline(xt,yt,xt1);
plot(xt1,yt1)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
title({'Joint 2 Angular Velocity';'(desired speed Vx = 2m/s, Vy = 0)'})

subplot(3,3,3)
xt = 1:20;
yt = phid(3,:);
xt1 = 1:0.1:20;
yt1 = spline(xt,yt,xt1);
plot(xt1,yt1)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
title({'Joint 3 Angular Velocity';'(desired speed Vx = 2m/s, Vy = 0)'})

subplot(3,3,4)
xt = 1:20;
yt = phid(4,:);
xt1 = 1:0.1:20;
yt1 = spline(xt,yt,xt1);
plot(xt1,yt1)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
title({'Joint 4 Angular Velocity';'(desired speed Vx = 2m/s, Vy = 0)'})

subplot(3,3,5)
xt = 1:20;
yt = phid(5,:);
xt1 = 1:0.1:20;
yt1 = spline(xt,yt,xt1);
plot(xt1,yt1)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
title({'Joint 5 Angular Velocity';'(desired speed Vx = 2m/s, Vy = 0)'})

subplot(3,3,6)
xt = 1:20;
yt = phid(6,:);
xt1 = 1:0.1:20;
yt1 = spline(xt,yt,xt1);
plot(xt1,yt1)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
title({'Joint 6 Angular Velocity';'(desired speed Vx = 2m/s, Vy = 0)'})

subplot(3,3,7)
xt = 1:20;
yt = phid(7,:);
xt1 = 1:0.1:20;
yt1 = spline(xt,yt,xt1);
plot(xt1,yt1)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
title({'Joint 7 Angular Velocity';'(desired speed Vx = 2m/s, Vy = 0)'})


