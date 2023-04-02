function animate_mov_data(data)

LimbNames={'Head','UpperarmR','ForearmR','HandR','Thorax','Abdomen','HipR','LegR','FootR','FootL','LegL','HipL','HandL','ForearmL','UpperarmL'};
JointNames={'HeadJoint','ShoulderJointR','ElbowJointR','WristJointR','ThoraxJoint','AbdomenJoint','HipJointR','KneeJointR','AnkleJointR','AnkleJointL','KneeJointL','HipJointL','WristJointL','ElbowJointL','ShoulderJointL'};
LimbIndex=[4,5,7,9,3,2,11,13,15,16,14,12,10,8,6];
limb_numbers=[1:15];

try
    close(W)
    delete(W)
catch
end
W=vrworld('W_trim_nofarch.wrl');
open(W);
Wfig=vrfigure(W)
dx=-1;dy=2.0;dz=1.5;

%W.Body2_pos.translation=[20 0 0];
W.Body_pos.translation=[0 0 0];
cam_pos=[dx, dy, dz];
W.current_viewpoint.position=cam_pos;

[m, n]=size(data);
time=0;
dt=1/240;
for j=1:1:n
    
    all_angs1=reshape(data(:,j),[3,15]);
    DOF=DOF_from_all_angs_Body(all_angs1);
    
    psi=DOF(LimbIndex(limb_numbers),1);
    theta=DOF(LimbIndex(limb_numbers),2);
    phi=DOF(LimbIndex(limb_numbers),3);
    for i=1:length(limb_numbers)
        [q, vec_euler, ang_euler]=calc_quaternion_from_euler_angles(['y','z','x'],[psi(i), -theta(i), phi(i)]);
        str=['W.' JointNames{limb_numbers(i)} '_pos.rotation=[vec_euler, ang_euler];'];
        eval(str);
    end
    
    drawnow
    %pause
    set(W,'Time',time);
    time=time+dt;
end
