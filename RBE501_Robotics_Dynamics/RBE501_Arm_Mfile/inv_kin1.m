function qs = inv_kin1(xd,yd,zd,theta)
l1 = 0.039;
l2 = 0.071;
l3 = 0.069;
l4 = 0.076;
le = 0;
c3 = sqrt(xd^2+yd^2+(zd-l1)^2);
if c3<(l2+l3+l4+le)
    q1 = atan2(yd,xd);
    q5 = 0;
    x3 = (sqrt(xd^2+yd^2)-cos(theta)*(l4+le))/sqrt(xd^2+yd^2)*xd;
    y3 = (sqrt(xd^2+yd^2)-cos(theta)*(l4+le))/sqrt(xd^2+yd^2)*yd;
    z3 = zd+sin(theta)*(l4+le);
    c2 = sqrt(x3^2+y3^2+(z3-l1)^2);
    c1 = sqrt(x3^2+y3^2+z3^2);
    q2 = pi-acos((l1^2+c2^2-c1^2)/(2*l1*c2))-acos((l2^2+c2^2-l3^2)/(2*l2*c2));
    q3 = 3*pi/4-acos((l2^2+l3^2-c2^2)/(2*l2*l3));
    q4 = 3*pi/4-acos((c2^2+(l4+le)^2-c3^2)/(2*c2*(l4+le)))-...
         acos((l3^2+c2^2-l2^2)/(2*l3*c2));
    qs = [q1;q2;q3;q4;q5];
end
