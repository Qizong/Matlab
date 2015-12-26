% -----------------------------------------------------
% -----------------------------------------------------
% Title:    cubic
%           [Compute the Cubic Polynomial trajectory for a joint variable]
% -----------------------------------------------------
% Inputs:   t0 = initial time
%           q0 = initial joint position
%           v0 = initial joint velocity
%           tf = final time
%           qf = final joint position
%           vf = final joint velocity
%
% Outputs:  X = Cubic Coefficients
% -----------------------------------------------------
% Useage: X=cubic(0,10,0,1,-20,0)
% -----------------------------------------------------
function X=cubic(t0,q0,v0,tf,qf,vf)

    %Write equations in matrix form
    A=[t0^3 t0^2 t0 1   
       3*t0^2 2*t0 1 0   
       tf^3 tf^2 tf 1  
       3*tf^2  2*tf 1 0];
    B=[q0;v0;qf;vf]; %column vector

    %Solve for polynomial coefficients
    %A*X=B
    X=A\B;

    %Define time vector
    t=linspace(t0,tf); %default is 100 points

    %Joint Space Trajectory
    q=X(1)+X(2)*t+X(3)*t.^2+X(4)*t.^3;
    qd=X(2)+2*X(3)*t+3*X(4)*t.^2;
    qdd=2*X(3)+6*X(4)*t;

 end
