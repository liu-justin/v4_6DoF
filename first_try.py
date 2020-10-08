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

# inverse dynamics
# need Glist, or spatial inertai matrix list
#   6x6 matrix, top left corner is 3x3 rotational inertia matrix, bottom right is mass of link * identity matrix
GList = []
np.set_printoptions(precision=7, suppress=True)
for link in obj.robot.link[2:-2]: # need to skip type fixed
    mass = float(link.inertial.mass["value"])

# got these values from solidworks, but have no idea why they are different than w,v
    # ix = [float(n) for n in link.inertial.origin["ix"].split()]
    # iy = [float(n) for n in link.inertial.origin["iy"].split()]
    # iz = [float(n) for n in link.inertial.origin["iz"].split()]
    # principle_axes = np.c_[ix,iy,iz]

    # translate from parent joint to CoM 
    # negative one b/c this is from parent link origin to CoM, but I need CoM to parent link origin
    xyz_CoM= -1*np.array([float(n) for n in link.inertial.origin["xyz"].split()])

    # grab Ixx, Ixy, Ixz, Iyy, Iyz, Izz about the CoM, with the parent link coordinate systems
    inertia_values_CoM = [float(n) for n in (vars(link.inertial.inertia_CoM)["_attributes"].values())]
    
    # putting those values into a rotational inertia matrix, centered at CoM, using parent link coords
    I_CoM = np.array([[inertia_values_CoM[0], inertia_values_CoM[1], inertia_values_CoM[2]],
                      [inertia_values_CoM[1], inertia_values_CoM[3], inertia_values_CoM[4]],
                      [inertia_values_CoM[2], inertia_values_CoM[4], inertia_values_CoM[5]]])

    # grabbing the eigenvectors of the rotational inertia matrix, to find the principle axes of inertia
    w,v = np.linalg.eig(I_CoM)
    print(f"eigenvectors:\n {v}")  
    print(f"inertia: \n{I_CoM}")
    
    # rotational inertia matrix, centered at CoM, aligned w/ principle axes of inertia
    rotated_I_CoM = np.transpose(v) @ I_CoM @ v
    print(f"inertia about rotated coords: \n{rotated_I_CoM}")
    # rotational inertia matrix, centered at parent link origin, aligned w/ parent link origin coords
    translated_T_CoM = I_CoM + mass*(np.inner(xyz_CoM, xyz_CoM)*np.identity(3) - np.outer(xyz_CoM, xyz_CoM))
    print(f"inertial rotational matrix at parent link: \n{translated_T_CoM}")

    # translated_T_CoM is pretty close to the value obtained from SOLIDWORKS
    # inertia_values_joint = [float(n) for n in (vars(link.inertial.inertia_joint)["_attributes"].values())]
    # I_joint = np.array([[inertia_values_joint[0], inertia_values_joint[1], inertia_values_joint[2]],
    #                     [inertia_values_joint[1], inertia_values_joint[3], inertia_values_joint[4]],
    #                     [inertia_values_joint[2], inertia_values_joint[4], inertia_values_joint[5]]])
 
    mI = mass*np.identity(3)
    zeros = np.zeros((3,3))
    Gi = np.c_[np.r_[rotated_I_CoM, zeros], np.r_[zeros,mI]]
    GList.append(Gi)

    print(f"Gi:\n{np.round(Gi,5)}")

print(TList)

# crete a trajectory to follow using functions from Ch9
theta_start = np.array([0,0,0,0,0,0])
theta_end = np.array([0,-1*np.pi/2, np.pi/2, 0,0,0])

# final time
T_final = 3
N = 1000
method = 5
# use a joint trajectory to get from rest to home
# then use a cartesian trajectory to get other places
# create the trajectory, N x n matrix where each row is  n-vector of joint variables at an instant in time
trajectory = mr.CartesianTrajectory(theta_start, theta_end, T_final, N, method)
theta_matrix = np.array(trajectory).copy()
d_theta_matrix = np.zeros((1000,3))
dd_theta_matrix = np.zeros((1000,3))
