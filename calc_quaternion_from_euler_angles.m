function [q, vec_euler, ang_euler]=calc_quaternion_from_euler_angles(seq_rotaxis,seq_rotvalues)

%any sequence of euler rotations: seq_rotaxis['x', 'y', 'z'] seq_values=[pi/4 pi/3 pi/6]

q_all=zeros(4,length(seq_rotvalues));
%calc individual rotation quaternions
for i=1:length(seq_rotvalues)
    switch seq_rotaxis(i)
        case 'x'
            q_all(:,i)=[cos(seq_rotvalues(i)/2);sin(seq_rotvalues(i)/2);0;0];
        case 'y'
            q_all(:,i)=[cos(seq_rotvalues(i)/2);0;sin(seq_rotvalues(i)/2);0];
        case 'z'
            q_all(:,i)=[cos(seq_rotvalues(i)/2);0;0;sin(seq_rotvalues(i)/2)];
        otherwise
    end
end


q=q_all(:,1);
for i=2:length(seq_rotvalues)
    q_mat=[q(1) -q(2) -q(3) -q(4); q(2) q(1) -q(4) q(3); q(3) q(4) q(1) -q(2); q(4) -q(3) q(2) q(1)];
    q=q_mat*q_all(:,i);
end

ang_euler=2*acos(q(1));
vec_euler=zeros(1,3);
if ang_euler~=0
vec_euler(1)=q(2) / sqrt(1-q(1)^2);
vec_euler(2)=q(3) / sqrt(1-q(1)^2);
vec_euler(3)=q(4) / sqrt(1-q(1)^2);
end