## Project: Kinematics Pick & Place
### Project 2 submission by Lewis Siempelkamp for Udacity's Robotics NanoDegree Term 1 
### February 2019
---


**Steps to complete the project:**  


1. Set up your ROS Workspace.
2. Download or clone the [project repository](https://github.com/udacity/RoboND-Kinematics-Project) into the ***src*** directory of your ROS Workspace.  
3. Experiment with the forward_kinematics environment and get familiar with the robot.
4. Launch in [demo mode](https://classroom.udacity.com/nanodegrees/nd209/parts/7b2fd2d7-e181-401e-977a-6158c77bf816/modules/8855de3f-2897-46c3-a805-628b5ecf045b/lessons/91d017b1-4493-4522-ad52-04a74a01094c/concepts/ae64bb91-e8c4-44c9-adbe-798e8f688193).
5. Perform Kinematic Analysis for the robot following the [project rubric](https://review.udacity.com/#!/rubrics/972/view).
6. Fill in the `IK_server.py` with your Inverse Kinematics code. 


[//]: # (Image References)

[image1]: ./misc_images/img_1.jpg
[image2]: ./misc_images/img_2.jpg
[image3]: ./misc_images/misc2.png

## [Rubric](https://review.udacity.com/#!/rubrics/972/view) Points
### Here I will consider the rubric points individually and describe how I addressed each point in my implementation.  

---
### Writeup / README

#### 1. Provide a Writeup / README that includes all the rubric points and how you addressed each one.  You can submit your writeup as markdown or pdf.  

You're reading it!

### Kinematic Analysis
#### 1. Run the forward_kinematics demo and evaluate the kr210.urdf.xacro file to perform kinematic analysis of Kuka KR210 robot and derive its DH parameters.

By parsing the supplied [URDF file](kuka_arm/urdf) for the Kuka Arm the joint coordinate frames can be extracted and translated into coordinate frames following the modified DH convention 

![DH Params Sketch][image1]

Links | alpha(i-1) | a(i-1) | d(i) | theta(i)
--- | --- | --- | --- | ---
0->1 | 0 | 0 | 0.75 | q1
1->2 | - pi/2 | 0.35 | 0 | -pi/2 + q2
2->3 | 0 | 1.25 | 0 | q3
3->4 |  -pi/2 | -0.054 | 1.5 | q4
4->5 | pi/2 | 0 | 0 | q5
5->6 | -pi/2 | 0 | 0 | q6
6->EE | 0 | 0 | 0.303 | 0

#### 2. Using the DH parameter table you derived earlier, create individual transformation matrices about each joint. In addition, also generate a generalized homogeneous transform between base_link and gripper_link using only end-effector(gripper) pose.

From the DH table a generalized Homogeneous transform from the base to the EE can be calculated using the DH method:

#Define Modified DH Transformation Matrix
        
        def TF_Matrix(alpha, a, d, q):
            TF = Matrix([[            cos(q),           -sin(q),           0,             a],
                        [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                        [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                        [                 0,                 0,           0,             1]])
            return TF
        
        #Create individual transformation matrices
        T0_1 =  TF_Matrix(alpha0, a0, d1, q1).subs(DH_Table)
        T1_2 =  TF_Matrix(alpha1, a1, d2, q2).subs(DH_Table)
        T2_3 =  TF_Matrix(alpha2, a2, d3, q3).subs(DH_Table)
        T3_4 =  TF_Matrix(alpha3, a3, d4, q4).subs(DH_Table)
        T4_5 =  TF_Matrix(alpha4, a4, d5, q5).subs(DH_Table)
        T5_6 =  TF_Matrix(alpha5, a5, d6, q6).subs(DH_Table)
        T6_EE = TF_Matrix(alpha6, a6, d7, q7).subs(DH_Table)

        T0_EE = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_EE


The DH conventional coordinate frames do not match with the coordinate frames supplied in the URDF, however, so the final transofrm must be rotated to align itself with the EE coord frame as defined in the URDF:

#Define RPY rotation matrices
    # http://planning.cs.uiuc.edu/node102.html

    r,p,y = symbols('r p y') #roll pitch yaw

    ROT_x = Matrix([[      1,      0,      0],
                    [      0, cos(r),-sin(r)],
                    [      0, sin(r), cos(r)]]) # ROLL
    
    ROT_y = Matrix([[ cos(p),      0, sin(p)],
                    [      0,      1,      0],
                    [-sin(p),      0, cos(p)]]) # PITCH

    ROT_z = Matrix([[ cos(y),-sin(y),      0],
                    [ sin(y), cos(y),      0],
                    [      0,      0,      1]]) # YAW

    ROT_EE = ROT_z * ROT_y * ROT_x
    
    # Generate rotation error correction for the EE between the URDF and DH params
    # by performing the following intrinsic rotation about EE z and y axes
    Rot_Error = ROT_z.subs(y,radians(180))*ROT_y.subs(p, radians(-90))
    ROT_EE = ROT_EE*Rot_Error
    T0_EE = T0_EE*(Rot_Error.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]])))
    
    
#### 3. Decouple Inverse Kinematics problem into Inverse Position Kinematics and inverse Orientation Kinematics; doing so derive the equations to calculate all individual joint angles.

To generate the Inverse Kinematics (IK) of this system. The joint angles (theta) must be extracted for any arbitrary End Effector (EE) Positiona and Orientation. This is done by decopuling the the Wrist Center (WC) Position from its orientation.
The WC position is extracted from the EE as follows: 

        WC = EEpos - d7*ROT_EE[:,2]# WC Position Matrix

With the WC position calculated, theta 1, theta2, and theta3 can be calculated using some trig:

![Sketch of WC Position Calcs][image2]

With thetas 1-3 calculated the orientation can be determined by realizing that the rotation of the EE is equal to the rotation of the the Wrist, that is that it is equal to the Rotational Transform from Joint 3 to joint 6 which can be calculated as follows:

#extract rotation matrix from base to WC

    R0_3 = T0_1[0:3,0:3] * T1_2[0:3,0:3] * T2_3[0:3,0:3]
    R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
    #extract rotation matrix from WC to EE
    R3_6 = R0_3.transpose() * ROT_EE #(where ROT_EE is known and has been corrected to the DH coordinate frame convention)

By subbing in the joint anngles theta1-3 into R0_3.transpose the RHS of the equation, R3_6, can be simplified and used to solve for the remaining joint angles, theta4-6:

#calculate the Wrist Joint angles Joints 3-6

    theta4 = atan2(R3_6[2,2], -R3_6[0,2])
    theta5 = atan2(sqrt(R3_6[0,2]**2 + R3_6[2,2]**2),R3_6[1,2])
    theta6 = atan2(-R3_6[1,1],R3_6[1,0])

Now we are able to solve for all of the joint angles given any arbitrary EE Position and Orientation and have solved the IK for this arm. 


### Project Implementation

#### 1. Fill in the `IK_server.py` file with properly commented python code for calculating Inverse Kinematics based on previously performed Kinematic Analysis. Your code must guide the robot to successfully complete 8/10 pick and place cycles. Briefly discuss the code you implemented and your results. 

The theory above was codified and debugged initially using [IK_Debug.py](/IK_Debug) to ensure the FK and IK were solved correctly.

By launching the provided RViz environment, ($roslaunch kuka_arm forward_kinematics.launch),  you can set the Kuka Arm to arbitrary joint angles and extract the EE position/orientation which can then be written into the IK_Debug.py file to verify your FK and IK results in the expected theta values.

Once certain the FK and IK work as expected they can be ported into the [IK_Server.py](/kuka_arm/scripts/)


And just for fun, another example image:
![alt text][image3]


