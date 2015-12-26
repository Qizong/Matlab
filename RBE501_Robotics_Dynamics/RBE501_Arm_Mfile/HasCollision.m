function R = HasCollision(q1, q2, q3, q4, q5)
% pass in the four arm joints, only check the first link and the third link
q1 = q1/pi*180;
q2 = q2/pi*180;
q3 = q3/pi*180;
q4 = q4/pi*180;
q5 = q5/pi*180;

l1 = 0.039;
l2 = 0.071;
l3 = 0.069;
l4 = 0.076;
le = 0;

J1x = 0;
J1y = l1;

J2x = J1x + l2*cosd(q2);
J2y = J1y + l2*sind(q2);

J3x = J2x + l3*cosd(q2+q3);
J3y = J2y + l3*sind(q2+q3);

J4x = J3x + l4*cosd(q2+q3+q4);
J4y = J3y + l4*sind(q2+q3+q4);

L1 = [J1x J2x; J1y J2y];
L2 = [J2x J3x; J2y J3y];
L3 = [J3x,J4x; J3y J4y];

% figure(1);
% hold on;
% axis equal;
% line(L1(1,:),L1(2,:));
% line(L2(1,:),L2(2,:));
% line(L3(1,:),L3(2,:));
% hold off;
% check if two link has collision
% Return 1 when intersect
% Return 0 when no intersection
if size(CheckIt(L1,L3),2) == 0
    R = 0;
else
    R = 1; 
end