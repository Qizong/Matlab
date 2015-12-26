function RBE501_arm_main()
clear all
clc
%% 1 - D-H Parameters
% Set symbolic parameters for robot
l1 = 0.039;
l2 = 0.071;
l3 = 0.069;
l4 = 0.076;
le = 0;
syms q1 q2 q3 q4 q5
% D-H table
% Notation 'alphai' is the rotation around the x-axis, 'ai' is the 
% translation along the x-axis, 'thetai' is the rotation around the z-axis,
% and 'di' is the translation along the z-axis.

ai = [0 l2 l3 0 0];
alphai = [-pi/2 0 0 -pi/2 0];
alphai_d = alphai/pi*180;
di = [l1 0 0 0 l4+le];
thetai = [q1 q2-pi/2 q3+pi/4 q4-pi/4 q5];

%% 2 - Arm Forward Position Kinematics
% Transformation matrices
T_01 = RZ(thetai(1)) * TZ(di(1)) * TX(ai(1)) * RdX(alphai_d(1))
T_12 = RZ(thetai(2)) * TZ(di(2)) * TX(ai(2)) * RdX(alphai_d(2))
T_23 = RZ(thetai(3)) * TZ(di(3)) * TX(ai(3)) * RdX(alphai_d(3))
T_34 = RZ(thetai(4)) * TZ(di(4)) * TX(ai(4)) * RdX(alphai_d(4))
T_4e = RZ(thetai(5)) * TZ(di(5)) * TX(ai(5)) * RdX(alphai_d(5))
% Transformation matrix from Frame R to Frame T
T_0e = T_01 * T_12 * T_23 * T_34 * T_4e
% Tip position
Px = T_0e(1,4)
Py = T_0e(2,4)
Pz = T_0e(3,4)

%% 3 - Initial Position Plot
q1_0 = 0;
q2_0 = 0;
q3_0 = 0;
q4_0 = 0;
q5_0 = 0;
figure(1)
T = plot_arm(q1_0,q2_0,q3_0,q4_0,q5_0);
title('Initial Position')
[pp1,pp2,pp3,pp4,ppe] = for_kin(0,0,0,0,0);
ppe

%% 4 - Inverse Position Kinematics Example
qs_1 = inv_kin(0.06,0.06,0.02)
figure(2)
T1 = plot_arm(qs_1(1),qs_1(2),qs_1(3),qs_1(4),qs_1(5));
title('Arm Configuration with Desired Positoin [0.06 0.06 0.02](m)')

%% 5 - Rectilinear Trajectory
% Desired tajectory: from [0.06,0.06,0.02] to [0.06,-0.06,0.02]
x0 = 0.12; y0 = 0.1; z0 = 0.1288;
xf = 0.12; yf = -0.1; zf = 0.1288;
n = 10;
for i = 1:n+1
    xtemp = x0;
    ytemp = y0-(i-1)/n*(y0-yf);
    ztemp = z0;
    qs(:,i) = inv_kin1(xtemp,ytemp,ztemp,0);
    R = HasCollision(qs(1,i),qs(2,i),qs(3,i),qs(4,i),qs(5,i));
    if R == 1
        fprintf('There is collision!')
        break
    end
end


qj1(:,1) = qs(1,:);
qj1(:,2) = linspace(0.1,0.1,n+1);
qj1(:,3) = linspace(0,0,n+1);
qj2(:,1) = qs(2,:);
qj2(:,2) = linspace(0.1,0.1,n+1);
qj2(:,3) = linspace(0,0,n+1);
qj3(:,1) = qs(3,:);
qj3(:,2) = linspace(0.1,0.1,n+1);
qj3(:,3) = linspace(0,0,n+1);
qj4(:,1) = qs(4,:);
qj4(:,2) = linspace(0.1,0.1,n+1);
qj4(:,3) = linspace(0,0,n+1);
qj5(:,1) = qs(5,:);
qj5(:,2) = linspace(0.1,0.1,n+1);
qj5(:,3) = linspace(0,0,n+1);

ts = 1; te = 3;
J1_in1 = timeseries(qj1,ts:(te-ts)/n:te);
J2_in1 = timeseries(qj2,ts:(te-ts)/n:te);
J3_in1 = timeseries(qj3,ts:(te-ts)/n:te);
J4_in1 = timeseries(qj4,ts:(te-ts)/n:te);
J5_in1 = timeseries(qj5,ts:(te-ts)/n:te);
save('AngleSet1_line','J1_in1','J2_in1','J3_in1','J4_in1','J5_in1');

J_line(:,1) = J1_in1.Time;
J_line(:,2) = J1_in1.Data(:,1)/pi*180;
J_line(:,3) = J2_in1.Data(:,1)/pi*180;
J_line(:,4) = J3_in1.Data(:,1)/pi*180+45;
J_line(:,5) = J4_in1.Data(:,1)/pi*180+45;
J_line(:,6) = J5_in1.Data(:,1)/pi*180;
J_line_xls = [0 0 0 0 0 0;J_line];

xlswrite('out_line.xls',J_line_xls)
% J1_line
%xlswrite('lineJ1.xlsx',J1_in1)

figure(3)
subplot(3,2,1)
plot(0:1/n:1,qs(1,:))
title('Joint 1 Angle')
xlabel('Time (s)')
ylabel('Joint Angle (rad)')

subplot(3,2,2)
plot(0:1/n:1,qs(2,:))
title('Joint 2 Angle')
xlabel('Time (s)')
ylabel('Joint Angle (rad)')

subplot(3,2,3)
plot(0:1/n:1,qs(3,:))
title('Joint 3 Angle')
xlabel('Time (s)')
ylabel('Joint Angle (rad)')

subplot(3,2,4)
plot(0:1/n:1,qs(4,:))
title('Joint 4 Angle')
xlabel('Time (s)')
ylabel('Joint Angle (rad)')

subplot(3,2,5)
plot(0:1/n:1,qs(5,:))
title('Joint 5 Angle')
xlabel('Time (s)')
ylabel('Joint Angle (rad)')


%% 6 - Coefficient of Cubic Polynomial Trajectory
t0 = 0; tf = 1; v0 = 0; vf = 0.01; q0 = 0;
qs_1 = inv_kin1(0.06,0.06,0.02,pi/3);
A1 = cubic(t0,q0,v0,tf,qs_1(1),vf);
A2 = cubic(t0,q0,v0,tf,qs_1(2),vf);
A3 = cubic(t0,q0,v0,tf,qs_1(3),vf);
A4 = cubic(t0,q0,v0,tf,qs_1(4),vf);
A5 = cubic(t0,q0,v0,tf,qs_1(5),vf);
save('Cubic','A1','A2','A3','A4','A5')

%% 7 - Cirle Trajectory

% Center of circle
x0 = 0.1488; y0 = 0; z0 = 0.1288;
r = 0.02; 
n = 50;
for i = 1:n+1
    theta_circle = (i-1)*2*pi/n;
    x_cl = x0;
    z_cl = z0+r*cos(theta_circle);
    y_cl = y0+r*sin(theta_circle);
    qs_cl(:,i) = inv_kin1(x_cl,y_cl,z_cl,0);
end

%time 1 to 5
qj1_cl(:,1) = qs_cl(1,:);
qj1_cl(:,2) = linspace(0.1,0.1,n+1);
qj1_cl(:,3) = linspace(0,0,n+1);
qj2_cl(:,1) = qs_cl(2,:);
qj2_cl(:,2) = linspace(0.1,0.1,n+1);
qj2_cl(:,3) = linspace(0,0,n+1);
qj3_cl(:,1) = qs_cl(3,:);
qj3_cl(:,2) = linspace(0.1,0.1,n+1);
qj3_cl(:,3) = linspace(0,0,n+1);
qj4_cl(:,1) = qs_cl(4,:);
qj4_cl(:,2) = linspace(0.1,0.1,n+1);
qj4_cl(:,3) = linspace(0,0,n+1);
qj5_cl(:,1) = qs_cl(5,:);
qj5_cl(:,2) = linspace(0.1,0.1,n+1);
qj5_cl(:,3) = linspace(0,0,n+1);qj1_cl(:,1) = qs_cl(1,:);
%time 0 to 1
n1 = 20;
qj1_cl0(:,1) = linspace(0,qs_cl(1,1),n1+1);
qj1_cl0(:,2) = linspace(0.1,0.1,n1+1);
qj1_cl0(:,3) = linspace(0,0,n1+1);
qj2_cl0(:,1) = linspace(0,qs_cl(2,1),n1+1);
qj2_cl0(:,2) = linspace(0.1,0.1,n1+1);
qj2_cl0(:,3) = linspace(0,0,n1+1);
qj3_cl0(:,1) = linspace(0,qs_cl(3,1),n1+1);
qj3_cl0(:,2) = linspace(0.1,0.1,n1+1);
qj3_cl0(:,3) = linspace(0,0,n1+1);
qj4_cl0(:,1) = linspace(0,qs_cl(4,1),n1+1);
qj4_cl0(:,2) = linspace(0.1,0.1,n1+1);
qj4_cl0(:,3) = linspace(0,0,n1+1);
qj5_cl0(:,1) = linspace(0,qs_cl(5,1),n1+1);
qj5_cl0(:,2) = linspace(0.1,0.1,n1+1);
qj5_cl0(:,3) = linspace(0,0,n1+1);
%generate time series (0 to 1)
n1 = 20;
ts = 0; te = 1;
J1_cl0 = timeseries([qj1_cl0],ts:(te-ts)/n1:te);
J2_cl0 = timeseries([qj2_cl0],ts:(te-ts)/n1:te);
J3_cl0 = timeseries([qj3_cl0],ts:(te-ts)/n1:te);
J4_cl0 = timeseries([qj4_cl0],ts:(te-ts)/n1:te);
J5_cl0 = timeseries([qj5_cl0],ts:(te-ts)/n1:te);
save('AngleSet1_circle(0to1s)','J1_cl0','J2_cl0','J3_cl0','J4_cl0','J5_cl0');
%generate time series (1 to 5)
ts = 1; te = 5;
J1_cl = timeseries([qj1_cl],ts:(te-ts)/n:te);
J2_cl = timeseries([qj2_cl],ts:(te-ts)/n:te);
J3_cl = timeseries([qj3_cl],ts:(te-ts)/n:te);
J4_cl = timeseries([qj4_cl],ts:(te-ts)/n:te);
J5_cl = timeseries([qj5_cl],ts:(te-ts)/n:te);
save('AngleSet1_circle1(1to5s)','J1_cl','J2_cl','J3_cl','J4_cl','J5_cl');

J_cl(:,1) = J1_cl.Time;
J_cl(:,2) = J1_cl.Data(:,1)*180/pi;
J_cl(:,3) = J2_cl.Data(:,1)*180/pi;
J_cl(:,4) = J3_cl.Data(:,1)*180/pi+45;
J_cl(:,5) = J4_cl.Data(:,1)*180/pi+45;
J_cl(:,6) = J5_cl.Data(:,1)*180/pi;

J_cl_xls = [0 0 0 0 0 0;J_cl];
xlswrite('out_circle.xls',J_cl_xls)


%% 8 - Board Trajectory (Star shape)
Inp1 = xlsread('out_test.xlsx');
%Inp1 = Inp1(1:length(Inp1(:,1))/10:length(Inp1(:,1)), :);

nn = length(Inp1(:,1))


for i = 1:nn
    qd(:,i) = inv_kin1(Inp1(i,1),Inp1(i,2),Inp1(i,3),0);
end


qd1(:,1) = qd(1,:);
qd1(:,2) = linspace(0.1,0.1,nn);
qd1(:,3) = linspace(0,0,nn);
qd2(:,1) = qd(2,:);
qd2(:,2) = linspace(0.1,0.1,nn);
qd2(:,3) = linspace(0,0,nn);
qd3(:,1) = qd(3,:)-pi/3;
qd3(:,2) = linspace(0.1,0.1,nn);
qd3(:,3) = linspace(0,0,nn);
qd4(:,1) = qd(4,:);
qd4(:,2) = linspace(0.1,0.1,nn);
qd4(:,3) = linspace(0,0,nn);
qd5(:,1) = qd(5,:);
qd5(:,2) = linspace(0.1,0.1,nn);
qd5(:,3) = linspace(0,0,nn);

% ts = 0; te = 1;
% J1b0 = timeseries([linspace(0,qd1(1,1),nn)],ts:(te-ts)/(nn-1):te);
% J2b0 = timeseries([linspace(0,qd2(1,1),nn)],ts:(te-ts)/(nn-1):te);
% J3b0 = timeseries([linspace(0,qd3(1,1),nn)],ts:(te-ts)/(nn-1):te);
% J4b0 = timeseries([linspace(0,qd4(1,1),nn)],ts:(te-ts)/(nn-1):te);
% J5b0 = timeseries([linspace(0,qd5(1,1),nn)],ts:(te-ts)/(nn-1):te);
% save('AngleSet1_board_0to1','J1b0','J2b0','J3b0','J4b0','J5b0');

ts = 1; te = 11;
J1b = timeseries([qd1],ts:(te-ts)/(nn-1):te);
J2b = timeseries([qd2],ts:(te-ts)/(nn-1):te);
J3b = timeseries([qd3],ts:(te-ts)/(nn-1):te);
J4b = timeseries([qd4],ts:(te-ts)/(nn-1):te);
J5b = timeseries([qd5],ts:(te-ts)/(nn-1):te);
save('AngleSet1_board_1to11','J1b','J2b','J3b','J4b','J5b');

J_b(:,1) = J1b.Time;
J_b(:,2) = J1b.Data(:,1)/pi*180;
J_b(:,3) = J2b.Data(:,1)/pi*180;
J_b(:,4) = J3b.Data(:,1)/pi*180+45;
J_b(:,5) = J4b.Data(:,1)/pi*180+45;
J_b(:,6) = J5b.Data(:,1)/pi*180;

J_b_xls = [0 0 0 0 0 0;J_b];
xlswrite('out_STAR.xls',J_b_xls)

%% Functions used in script
function Rx = RX(theta)
Rx = [1 0 0 0; 
      0 cos(theta) -sin(theta) 0;
      0 sin(theta) cos(theta) 0; 
      0 0 0 1];
  
function Rdx = RdX(theta)
Rdx = [1 0 0 0; 
      0 cosd(theta) -sind(theta) 0;
      0 sind(theta) cosd(theta) 0; 
      0 0 0 1];

function Ry = RY(theta)
Ry = [cos(theta) 0 sin(theta) 0;
      0 1 0 0;
      -sin(theta) 0 cos(theta) 0;
      0 0 0 1];

function Rz = RZ(theta)
Rz = [cos(theta) -sin(theta) 0 0; 
      sin(theta) cos(theta) 0 0; 
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