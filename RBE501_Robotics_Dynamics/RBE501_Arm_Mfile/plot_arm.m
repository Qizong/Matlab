% -----------------------------------------------------
% -----------------------------------------------------
% Date:     02/20/2015
% -----------------------------------------------------
% Title:    Plotarm Function
% -----------------------------------------------------
% Inputs:   q1 q2 q3 q4 q5 
% Outputs:  Transformation Matrix
% -----------------------------------------------------
% -----------------------------------------------------
function T = plot_arm(q1, q2, q3, q4, q5)
l1 = 0.039;
l2 = 0.071;
l3 = 0.069;
l4 = 0.076;
le = 0;

ai = [0 l2 l3 0 0];
alphai = [-pi/2 0 0 -pi/2 0]/pi*180;
di = [l1 0 0 0 l4+le];
thetai = [q1 q2-pi/2 q3+pi/4 q4-pi/4 q5]/pi*180;

T0e = eye(4);
for i = 1:5
    TT{i} = RZ(thetai(i)) * TZ(di(i)) * TX(ai(i)) * RX(alphai(i)); 
    T0e = T0e * TT{i};
    P(1:3,i) = T0e(1:3,4);
end

T = T0e;
x = [0 P(1,:)];
y = [0 P(2,:)];
z = [0 P(3,:)];

plot3(x,y,z,'LineWidth',3)
hold on
plot3(x,y,z,'ro','MarkerFaceColor','r')
hold off
axis([-0.2,0.2,-0.2,0.2,0,0.2])
title('Arm Configuration')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on

function Rx = RX(theta)
Rx = [1 0 0 0; 
      0 cosd(theta) -sind(theta) 0;
      0 sind(theta) cosd(theta) 0; 
      0 0 0 1];

function Ry = RY(theta)
Ry = [cosd(theta) 0 sind(theta) 0;
      0 1 0 0;
      -sind(theta) 0 cosd(theta) 0;
      0 0 0 1];

function Rz = RZ(theta)
Rz = [cosd(theta) -sind(theta) 0 0; 
      sind(theta) cosd(theta) 0 0; 
      0 0 1 0; 
      0 0 0 1];

function Tx = TX(d)
Tx = [1 0 0 d; 
      0 1 0 0; 
      0 0 1 0; 
      0 0 0 1];

function Ty = TY(d)
Ty = [1 0 0 0; 
      0 1 0 d; 
      0 0 1 0; 
      0 0 0 1];

function Tz = TZ(d)
Tz = [1 0 0 0; 
      0 1 0 0; 
      0 0 1 d; 
      0 0 0 1];