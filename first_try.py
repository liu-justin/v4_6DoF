import mr 
import untangle
import numpy as np

# first untagnle the xml file, then create the inputs for ch4 forward kinematics

obj = untangle.parse("6DoF_URDF.xml")
# [child["name"] for child in o.root.child]
TList = []
bodyList = []
print(type(obj.robot.joint))

distal_omega = np.array([0,0,0])
distal_p = np.array([0,0,0])

for joint in reversed(obj.robot.joint): # need to skip type fixed
    rpy = [float(n) for n in joint.origin["rpy"].split()]
    R = mr.RollPitchYawToRot(rpy[0], rpy[1], rpy[2])

    p = [float(n) for n in joint.origin["xyz"].split()]

    T = np.round(mr.RpToTrans(R,p),10)
    TList.insert(0,T)

    current_omega = [float(n) for n in joint.axis["xyz"].split()]
    current_v = np.dot(current_omega, distal_p)

    body_axis = np.r_[current_omega, current_v]
    bodyList.insert(0,body_axis)

    distal_p = -1* p

    


# alright, DH parameters use frames, and PoE only needs home and EE frames
# however, URDF has all frames b/c it needs them for mass crap in ch 8
# for FK, can just use home frame, ee frame, and calculate the body axes from the frames
# mr.FKinBody()