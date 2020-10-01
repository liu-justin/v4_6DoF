import mr 
import untangle
import numpy as np

# first untagnle the xml file, then create the inputs for ch4 forward kinematics

obj = untangle.parse("6DoF_URDF.xml")
# [child["name"] for child in o.root.child]
TList = []
bodyList = []

# grabbing the last joint (ee_joint) and its xyz
rpy_ee = [float(n) for n in obj.robot.joint[len(obj.robot.joint)-1].origin["rpy"].split()]
R_ee = mr.RollPitchYawToRot(rpy_ee[0],rpy_ee[1],rpy_ee[2])
p_ee = [float(n) for n in obj.robot.joint[len(obj.robot.joint)-1].origin["xyz"].split()]
T_ee = mr.RpToTrans(R_ee, p_ee)
print (T_ee)

joint_list = [joint for joint in obj.robot.joint if joint["type"]!="fixed"]

for joint in reversed(joint_list): # need to skip type fixed
    rpy = [float(n) for n in joint.origin["rpy"].split()]
    R = mr.RollPitchYawToRot(rpy[0], rpy[1], rpy[2])

    p = np.array([float(n) for n in joint.origin["xyz"].split()])

    # this T takes previous joint to current joint, or is current joint relative to prev joint
    # T_56, T_lower_higher
    T = np.round(mr.RpToTrans(R,p),10)
    TList.insert(0,T)
    # print(np.round(T,5))

    current_omega = [float(n) for n in joint.axis["xyz"].split()]
    current_omega_skewed = mr.VecToso3(current_omega)
    # (R_ee, p_ee) = mr.TransToRp(mr.TransInv(T_ee))
    (R_ee, p_ee) = mr.TransToRp(T_ee)
    current_v = np.dot(current_omega_skewed, p_ee)

    body_axis = np.r_[current_omega, current_v]
    print(body_axis)
    bodyList.insert(0,body_axis)

    T_ee = np.round(np.dot(T, T_ee),5)
    # print(T_ee)
    


# alright, DH parameters use frames, and PoE only needs home and EE frames
# however, URDF has all frames b/c it needs them for mass crap in ch 8
# for FK, can just use home frame, ee frame, and calculate the body axes from the frames
# mr.FKinBody()