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

joint_list = [joint for joint in obj.robot.joint if joint["type"]!="fixed"]

for joint in reversed(joint_list): # need to skip type fixed
    rpy = [float(n) for n in joint.origin["rpy"].split()]
    R = mr.RollPitchYawToRot(rpy[0], rpy[1], rpy[2])

    p = np.array([float(n) for n in joint.origin["xyz"].split()])

    # this T takes previous joint to current joint, or is current joint relative to prev joint
    # T_56, T_lower_higher
    T = mr.RpToTrans(R,p)
    TList.insert(0,T)
    # print(np.round(T,5))

    current_omega = [float(n) for n in joint.axis["xyz"].split()]
    current_omega_skewed = mr.VecToso3(current_omega)
    (R_ee, p_ee) = mr.TransToRp(mr.TransInv(T_ee))
    # (R_ee, p_ee) = mr.TransToRp(T_ee)
    current_v = -1*np.dot(current_omega_skewed, p_ee)

    body_axis = np.r_[current_omega, current_v]
    # print(body_axis)
    bodyList.insert(0,body_axis)

    T_ee = np.dot(T, T_ee)

# alright, DH parameters use frames, and PoE only needs home and EE frames
# however, URDF has all frames b/c it needs them for mass crap in ch 8
# for FK, can just use home frame, ee frame, and calculate the body axes from the frames
# angleList = np.array([0,-1*np.pi/2,np.pi/2,0,-1*np.pi/2,0])
angleList = np.array([0,-1*np.pi/2,np.pi/2,0,0,0])
bodyList = np.array(bodyList).T

final = mr.FKinBody(T_ee, bodyList, angleList)
# print(np.round(T_ee, 5))
print(np.round(final,5))

M = np.array([[-1, 0,  0, 0],
                      [ 0, 1,  0, 6],
                      [ 0, 0, -1, 2],
                      [ 0, 0,  0, 1]])
Blist = np.array([[0, 0, -1, 2, 0,   0],
                          [0, 0,  0, 0, 1,   0],
                          [0, 0,  1, 0, 0, 0.1]]).T
thetalist = np.array([np.pi / 2.0, 3, np.pi])

# inverse dynamics
# need Glist, or spatial inertai matrix list
#   6x6 matrix, top left corner is 3x3 rotational inertia matrix, bottom right is mass of link * identity matrix
Glist = []

for link in obj.robot.link[2:-2]: # need to skip type fixed
    mass = float(link.inertial.mass["value"])

    # ix = [float(n) for n in link.inertial.origin["ix"].split()]
    # iy = [float(n) for n in link.inertial.origin["iy"].split()]
    # iz = [float(n) for n in link.inertial.origin["iz"].split()]
    # principle_axes = np.c_[ix,iy,iz]

    CoM_xyz = [float(n) for n in link.inertial.origin["xyz"].split()]

    inertia = [float(n) for n in (vars(link.inertial.inertia)["_attributes"].values())]
    
    Ib = np.array([[inertia[0], inertia[1], inertia[2]],
                   [inertia[1], inertia[3], inertia[4]],
                   [inertia[2], inertia[4], inertia[5]]])

    print(f"eigenvalues and vectors:\n {np.linalg.eig(Ib)}")

    w,v = np.linalg.eig(Ib)

    np.set_printoptions(precision=7, suppress=True)
    print(f"inertia: \n{Ib}")
    print(f"inertia about rotated coords: \n{np.transpose(v) @ Ib @ v}")
    
    transformed_Ib = np.transpose(v) @ Ib @ v

    mI = mass*np.identity(3)

    zeros = np.zeros((3,3))
    Gi = np.c_[np.r_[transformed_Ib, zeros], np.r_[zeros,mI]]
    Glist.append(Gi)

    print(f"Gi:\n{np.round(Gi,5)}")


# says each reference frame {i} is attached to the CoM of each link i,i=1...n
# base frame is denoted {0}
# frame at end effector denoted {n+1} (i think this is link6)

# two ways to do things
# use the link CoM, and convert screw axes, forces from referenced joint to referenced link
# or, use the joint CoM, and use the bottom inertial values in the URDF

# i will probably try both to see if i get the same result