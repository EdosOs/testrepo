function animate_mov_component(neutral_pose, PC, recording)

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

if recording
    set(Wfig,'Record2D','on');
    set(Wfig,'Record2DCompressMethod','none');
    set(W,'Recording','on');
end
time=0;
time_f=10;
dt=0.01;
for j=1:1:time_f/dt
    contol_input=1.0*sin(2*pi*0.25*time);
    comp_sum=contol_input*PC+neutral_pose;
    all_angs1=reshape(comp_sum,3,15);
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

if recording
    set(W,'Recording','off');
end