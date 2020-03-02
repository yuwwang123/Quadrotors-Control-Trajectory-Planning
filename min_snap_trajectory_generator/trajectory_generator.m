% Min Jerk
function [desired_state] = trajectory_generator(t, qn, varargin)
% TRAJECTORY_GENERATOR: Turn a Dijkstra or A* path into a trajectory
% t: A scalar, specifying inquiry time
%
% varargin: variable number of input arguments. In the framework,
% this function will first (and only once!) be called like this:
%
% trajectory_generator([],[], 0, path)
%
% i.e. map = varargin{1} and path = varargin{2}.
%
% path: A N x 3 matrix where each row is (x, y, z) coordinate of a
% point in the path. N is the total number of points in the path
%
% This is when you compute and store the trajectory.
%
% Later it will be called with only t and qn as an argument, at
% which point you generate the desired state for point t.
%
persistent path
persistent Ts
persistent T_tot
persistent x_coeffs
persistent y_coeffs
persistent z_coeffs
if ~isempty(varargin)  
    path = varargin{2};
    T_tot =  1.3 * getPathLength(path);
    Ts = gen_Ts(path,T_tot);
    Ts(1)= -0.1;
    %x(t)
    x_coeffs = genQuintic(Ts, path, 1);
    y_coeffs = genQuintic(Ts, path, 2);
    z_coeffs = genQuintic(Ts, path, 3);    
    desired_state.pos = [0;0;0];
    desired_state.vel = [0;0;0];
    desired_state.acc = [0;0;0];
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
    return
end

if t >= Ts(end)
    pos = [path(end,1);path(end,2);path(end,3)];
    vel = [0;0;0];
    acc = [0;0;0];
    
else   
    %determin which interval t falls in
    m = 0;
    while t > Ts(m+1) 
       m = m+1;
    end
    m = max(1,m);
    m = min(length(Ts)-1, m);
    
    pos_x = dot(x_coeffs(m,:),[1  t  t^2  t^3   t^4  t^5]);
    vel_x = dot(x_coeffs(m,:),[0  1  2*t  3*t^2  4*t^3  5*t^4]);
    acc_x = dot(x_coeffs(m,:),[0  0  2     6*t  12*t^2  20*t^3]);

    pos_y = dot(y_coeffs(m,:),[1  t  t^2  t^3   t^4  t^5]);
    vel_y = dot(y_coeffs(m,:),[0  1  2*t  3*t^2  4*t^3  5*t^4]);
    acc_y = dot(y_coeffs(m,:),[0  0  2     6*t  12*t^2  20*t^3]);

    pos_z = dot(z_coeffs(m,:),[1  t  t^2  t^3   t^4  t^5]);
    vel_z = dot(z_coeffs(m,:),[0  1  2*t  3*t^2  4*t^3  5*t^4]);
    acc_z = dot(z_coeffs(m,:),[0  0  2     6*t  12*t^2  20*t^3]);

    pos = [pos_x;pos_y;pos_z];
    vel = [vel_x;vel_y;vel_z];
    acc = [acc_x;acc_y;acc_z];

    % use the "persistent" keyword to keep your trajectory around
    % inbetween function calls
end
yaw = 0;
yawdot = 0;

%
% When called without varargin (isempty(varargin) == true), compute
% and return the desired state here.
%
desired_state.pos = pos(:);
desired_state.vel = vel(:);
desired_state.acc = acc(:);
desired_state.yaw = yaw;
desired_state.yawdot = yawdot;

end


function coeffs = genQuintic(Ts, path, axis)
    m = size(path,1) - 1;  % total num of intervals along the path
    M = zeros(6*m, 6*m);
    RHS = zeros(6*m,1);
    % vel intermediate constraints row(1 ~ m-1)
    for i = 1:m-1
        M(i, 6*i-5:6*i) = vel_term(Ts(i+1));
        M(i, 6*i+1:6*i+6) = - vel_term(Ts(i+1));
    end
    % acc intermediate constraints row(m ~ 2m-2)
    for i = 1:m-1
        M(i+m-1, 6*i-5:6*i) = acc_term(Ts(i+1));
        M(i+m-1, 6*i+1:6*i+6) = - acc_term(Ts(i+1));
    end
    % accd intermediate constraints row(2m-1 ~ 3m-3)
    for i = 1:m-1
        M(i+2*m-2, 6*i-5:6*i) = accd_term(Ts(i+1));
        M(i+2*m-2, 6*i+1:6*i+6) = - accd_term(Ts(i+1));
    end
    
    % accdd intermediate constraints row(3m-2 ~ 4m-4)
    for i = 1:m-1
        M(i+3*m-3, 6*i-5:6*i) = accdd_term(Ts(i+1));
        M(i+3*m-3, 6*i+1:6*i+6) = - accdd_term(Ts(i+1));
    end
    
    % pos constraints row(4m-3 ~ 6m-4) 2m
    for i = 1:m
        M(2*i-1+ 4*m-4, 6*i-5:6*i) = pos_term(Ts(i));
        M(2*i+ 4*m-4,  6*i-5:6*i) =  pos_term(Ts(i+1));
        RHS(2*i-1+ 4*m-4) = path(i,axis);
        RHS(2*i+ 4*m-4) = path(i+1,axis);
    end
      
    M(6*m-3, 1:6) = vel_term(Ts(1));
    M(6*m-2, 1:6) = acc_term(Ts(1)); 
    M(6*m-1, 6*m-5:6*m) = vel_term(Ts(m+1)); 
    M(6*m,   6*m-5:6*m) = acc_term(Ts(m+1)); 
    
    RHS(6*m-3) = 0;
    RHS(6*m-2) = 0; 
    RHS(6*m-1) = 0; 
    RHS(6*m) = 0;
    
    coeffs_v = transpose(M\RHS);
    coeffs = zeros(m,6);
    for i = 1:m
        coeffs(i,:) = coeffs_v(6*i-5:6*i);
    end
end

function Ts = gen_Ts(path, T_tot)
    Ts = zeros(1,size(path,1));
    dist_tot =0;
    for i = 1:size(path,1)-1
        dist = norm(path(i,:)-path(i+1,:));
        Ts(i+1) = Ts(i) + dist;
        dist_tot = dist_tot + dist;
    end
    Ts = Ts * T_tot / dist_tot;    
end

function result = pos_term(t)
    result = [1  t  t^2  t^3   t^4  t^5];
end

function result = vel_term(t)
    result = [0  1  2*t  3*t^2  4*t^3  5*t^4];
end

function result = acc_term(t)
    result = [0  0  2     6*t  12*t^2  20*t^3];
end

function result = accd_term(t)
    result = [0  0  0  6   24*t  60*t^2];
end

function result = accdd_term(t)
    result = [0  0  0  0   24  120*t];
end

function length = getPathLength(path)
    length = 0;
    for i=1:size(path,1)-1
        length = length + norm(path(i,:)-path(i+1,:))
    end
end

% % Min Snap
% function [desired_state] = trajectory_generator(t, qn, varargin)
% % TRAJECTORY_GENERATOR: Turn a Dijkstra or A* path into a trajectory
% % t: A scalar, specifying inquiry time
% %
% % varargin: variable number of input arguments. In the framework,
% % this function will first (and only once!) be called like this:
% %
% % trajectory_generator([],[], 0, path)
% %
% % i.e. map = varargin{1} and path = varargin{2}.
% %
% % path: A N x 3 matrix where each row is (x, y, z) coordinate of a
% % point in the path. N is the total number of points in the path
% %
% % This is when you compute and store the trajectory.
% %
% % Later it will be called with only t and qn as an argument, at
% % which point you generate the desired state for point t.
% %
% 
% persistent path
% persistent Ts
% persistent T_tot
% persistent x_coeffs
% persistent y_coeffs
% persistent z_coeffs
% if isempty(t)  
%     path = varargin{2};
% %    path(end-1,:)=[];
% %   path(end-2,:)=[];
%     T_tot =  1.3 * getPathLength(path);
%     Ts = gen_Ts(path,T_tot);
%     Ts(1)= -0.1;
%     %x(t)
%     x_coeffs = gen_min_snap(Ts, path, 1);
%     y_coeffs = gen_min_snap(Ts, path, 2);
%     z_coeffs = gen_min_snap(Ts, path, 3);
%     desired_state.pos = [0;0;0];
%     desired_state.vel = [0;0;0];
%     desired_state.acc = [0;0;0];
%     desired_state.yaw = 0;
%     desired_state.yawdot = 0;
%     return
%     
% end
% 
% if t >= Ts(end)
%     pos = [path(end,1);path(end,2);path(end,3)];
%     vel = [0;0;0];
%     acc = [0;0;0];
%     
% else   
%     %determin which interval t falls in
%     m = 0;
%     while t > Ts(m+1) 
%        m = m+1;
%     end
%     m = max(1,m);
%     m = min(length(Ts)-1, m);
%     
%     pos_x = dot(x_coeffs(m,:),pos_term(t));
%     vel_x = dot(x_coeffs(m,:),vel_term(t));
%     acc_x = dot(x_coeffs(m,:),acc_term(t));
% 
%     pos_y = dot(y_coeffs(m,:),pos_term(t));
%     vel_y = dot(y_coeffs(m,:),vel_term(t));
%     acc_y = dot(y_coeffs(m,:),acc_term(t));
% 
%     pos_z = dot(z_coeffs(m,:),pos_term(t));
%     vel_z = dot(z_coeffs(m,:),vel_term(t));
%     acc_z = dot(z_coeffs(m,:),acc_term(t));
% 
%     pos = [pos_x;pos_y;pos_z];
%     vel = [vel_x;vel_y;vel_z];
%     acc = [acc_x;acc_y;acc_z];
% 
%     % use the "persistent" keyword to keep your trajectory around
%     % inbetween function calls
% end
% yaw = 0;
% yawdot = 0;
% 
% %
% % When called without varargin (isempty(varargin) == true), compute
% % and return the desired state here.
% %
% desired_state.pos = pos(:);
% desired_state.vel = vel(:);
% desired_state.acc = acc(:);
% desired_state.yaw = yaw;
% desired_state.yawdot = yawdot;
% 
% end
% 
% function coeffs = gen_min_snap(Ts, path, axis)
%     m = size(path,1) - 1;  % total num of intervals along the path
%     M = zeros(8*m, 8*m);
%     RHS = zeros(8*m,1);
%     % vel intermediate constraints row(1 ~ m-1)
%     for i = 1:m-1
%         M(i, 8*i-7:8*i) = vel_term(Ts(i+1));
%         M(i, 8*i+1:8*i+8) = - vel_term(Ts(i+1));
%     end
%     % acc intermediate constraints row(m ~ 2m-2)
%     for i = 1:m-1
%         M(i+m-1, 8*i-7:8*i) = acc_term(Ts(i+1));
%         M(i+m-1, 8*i+1:8*i+8) = - acc_term(Ts(i+1));
%     end
%     % accd intermediate constraints row(2m-1 ~ 3m-3)
%     for i = 1:m-1
%         M(i+2*m-2, 8*i-7:8*i) = accd_term(Ts(i+1));
%         M(i+2*m-2, 8*i+1:8*i+8) = - accd_term(Ts(i+1));
%     end
%     
%     % accdd intermediate constraints row(3m-2 ~ 4m-4)
%     for i = 1:m-1
%         M(i+3*m-3, 8*i-7:8*i) = accdd_term(Ts(i+1));
%         M(i+3*m-3, 8*i+1:8*i+8) = - accdd_term(Ts(i+1));
%     end
%     
%     % accddd intermediate constraints row(4m-3 ~ 5m-5)
%     for i = 1:m-1
%         M(i+4*m-4, 8*i-7:8*i) = accddd_term(Ts(i+1));
%         M(i+4*m-4, 8*i+1:8*i+8) = - accddd_term(Ts(i+1));
%     end
%     
%     % accdddd intermediate constraints row(5m-4 ~ 6m-6)
%     for i = 1:m-1
%         M(i+5*m-5, 8*i-7:8*i) = accdddd_term(Ts(i+1));
%         M(i+5*m-5, 8*i+1:8*i+8) = - accdddd_term(Ts(i+1));
%     end
%     
%     % pos constraints row(6m-5 ~ 8m-6) 2m
%     for i = 1:m
%         M(2*i-1+ 6*m-6, 8*i-7:8*i) = pos_term(Ts(i));
%         M(2*i+ 6*m-6,  8*i-7:8*i) =  pos_term(Ts(i+1));
%         RHS(2*i-1+ 6*m-6) = path(i,axis);
%         RHS(2*i+ 6*m-6) = path(i+1,axis);
%     end
%       
%     M(8*m-5, 1:8) = vel_term(Ts(1));
%     M(8*m-4, 1:8) = acc_term(Ts(1)); 
%     M(8*m-3, 1:8) = accd_term(Ts(1)); 
%     M(8*m-2, 8*m-7:8*m) = vel_term(Ts(m+1)); 
%     M(8*m-1, 8*m-7:8*m) = acc_term(Ts(m+1));
%     M(8*m,   8*m-7:8*m) = accd_term(Ts(m+1));
%     
%     RHS(8*m-5) = 0;
%     RHS(8*m-4) = 0;
%     RHS(8*m-3) = 0;
%     RHS(8*m-2) = 0; 
%     RHS(8*m-1) = 0; 
%     RHS(8*m) = 0;
%     
%     coeffs_v = transpose(M\RHS);
%     coeffs = zeros(m,8);
%     for i = 1:m
%         coeffs(i,:) = coeffs_v(8*i-7:8*i);
%     end
% end
% 
% function Ts = gen_Ts(path, T_tot)
%     Ts = zeros(1,size(path,1));
%     dist_tot =0;
%     for i = 1:size(path,1)-1
%         dist = norm(path(i,:)-path(i+1,:));
%         Ts(i+1) = Ts(i) + dist;
%         dist_tot = dist_tot + dist;
%     end
%     Ts = Ts * T_tot / dist_tot;    
% end
% 
% function result = pos_term(t)
%     result = [1  t  t^2  t^3   t^4  t^5  t^6  t^7];
% end
% 
% function result = vel_term(t)
%     result = [0  1  2*t  3*t^2  4*t^3  5*t^4  6*t^5  7*t^6];
% end
% 
% function result = acc_term(t)
%     result = [0  0  2     6*t  12*t^2  20*t^3  30*t^4  42*t^5];
% end
% 
% function result = accd_term(t)
%     result = [0  0  0  6   24*t  60*t^2  120*t^3  210*t^4];
% end
% 
% function result = accdd_term(t)
%     result = [0  0  0  0   24  120*t  360*t^2  840*t^3];
% end
% 
% function result = accddd_term(t)
%     result = [0  0  0  0   0  120  720*t  2520*t^2];
% end
% 
% function result = accdddd_term(t)
%     result = [0  0  0  0   0  0  720  5040*t];
% end
% 
% function length = getPathLength(path)
%     length = 0;
%     for i=1:size(path,1)-1
%         length = length + norm(path(i,:)-path(i+1,:))
%     end
% end
