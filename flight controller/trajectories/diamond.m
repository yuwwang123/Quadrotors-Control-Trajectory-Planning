function [desired_state] = diamond(t, qn)
% DIAMOND trajectory generator for a diamond

% =================== Your code goes here ===================
% You have to set the pos, vel, acc, yaw and yawdot variables
% NOTE: the simulator will spawn the robot to be at the
%       position you return for t == 0


    pos0 = [0;0;0];
    posf = [1/4;sqrt(2);sqrt(2)];
    [out, outd, outdd] = genQuintic(t, 0, 3, 0, norm(posf - pos0));


if (3<= t) && (t< 5)
    pos0 = [1/4; sqrt(2); sqrt(2)];
    posf = [1/2; 0; 2*sqrt(2)];
   [out, outd, outdd] = genQuintic(t, 3, 5, 0, norm(posf - pos0));
end  

if (5<= t) && (t< 7)
    pos0 = [1/2; 0; 2*sqrt(2)];
    posf = [3/4; -sqrt(2); sqrt(2)];
   [out, outd, outdd] = genQuintic(t, 5, 7, 0, norm(posf - pos0));
end  

if (7<= t) && (t< 8.5)
    pos0 = [3/4; -sqrt(2); sqrt(2)];
    posf = [1; 0; 0];
   [out, outd, outdd] = genQuintic(t, 7, 8.5, 0, norm(posf - pos0));
end  

if (t>8.5)
    pos0 = [1; 0; 0];
    posf = [0;0;0];
    out = 0; outd=0; outdd=0;
end
    
pos = pos0 + out * (posf-pos0)/norm(posf-pos0);    
vel = outd *(posf-pos0)/norm(posf-pos0);
acc = outdd *(posf-pos0)/norm(posf-pos0);

yaw = 0;
yawdot = 0;
 

% =================== Your code ends here ===================

desired_state.pos = pos(:);
desired_state.vel = vel(:);
desired_state.acc = acc(:);
desired_state.yaw = yaw;
desired_state.yawdot = yawdot;

end

function [out, outd, outdd] = genQuintic(t, t0, tf,  p0,v0,a0,pf,vf,af)

    M =[1    t0  t0^2    t0^3    t0^4    t0^5;
        0    1   2*t0    3*t0^2  4*t0^3  5*t0^4;
        0    0   2       6*t0    12*t0^2 20*t0^3;
        1    tf  tf^2    tf^3    tf^4    tf^5;
        0    1   2*tf    3*tf^2  4*tf^3  5*tf^4;
        0    0   2       6*tf    12*tf^2 20*tf^3];

    constraints = [p0; v0 ;a0; pf; vf ; af];  
    a = M\constraints;
    out     = a(1) + a(2)*t + a(3)*t^2 + a(4)*t^3 + a(5)*t^4 + a(6)*t^5;
    outd    = a(2) + 2*a(3)*t + 3*a(4)*t^2 + 4*a(5)*t^3 + 5*a(6)*t^4;
    outdd   = 2*a(3) + 6*a(4)*t + 12*a(5)*t^2 + 20*a(6)*t^3;
end
