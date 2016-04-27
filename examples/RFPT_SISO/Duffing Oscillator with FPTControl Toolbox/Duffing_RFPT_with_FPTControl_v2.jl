#Duffig Oscillator
#Load the toolbox
using FPTControl

# if  Adaptive parameter value is 1:it will be use Fixed Point Transformation Method, else (example 0) it will be just a simple PID Controller
Adaptive=1
#Simulation Lenght
LONG=Int(1e4)
# time step
δt=1e-3

#PID tune parameter
Λ=3

#Nominal Trajectory for a simple sin() type Nominal Trajectory
ω=5
Ampl=0.4

#Controller Parameters
K=1e5
B=-1
A=1e-5

# Exact model parameters
αₑ=1
δₑ=0.2
βₑ=1

# Approx. model parameters
αₐ=0.8
δₐ=0.1
βₐ=0.9

###############
# Exact Model #
###############

function exact(q,q_p,u)
  out=αₑ*q+δₑ*q_p+βₑ*q^3+u
  return out #q_pp
end

##########################
# Approx Model (inverse model)    #
##########################

function approx(q,q_p,q_pp)
  out=q_pp-αₐ*q-δₐ*q_p-βₐ*q^3
  return out #u
end

######################################
# Arrays for Plotting and simulation #
######################################

qN=zeros(LONG)
qN_p=zeros(LONG)
qN_pp=zeros(LONG)

qDes_pp=zeros(LONG)

qDef_pp=zeros(LONG)
u=zeros(LONG)

q=zeros(LONG)
q_p=zeros(LONG)
q_pp=zeros(LONG)

t=zeros(LONG)

#Initial condtions
error_int=0
past_input=0
past_response=0

#simulation

for i=1:LONG-1
 #store time
 t[i]=δt*i
 #Nominal Trajectory and the derivatives
 qN[i]=Ampl*sin(ω*t[i])
 qN_p[i]=Ampl*ω*cos(ω*t[i])
 qN_pp[i]=-Ampl*ω^2*sin(ω*t[i])
 #error (by definition)
 h=qN[i]-q[i]
 #first time derivative of the error
 h_p=qN_p[i]-q_p[i]
 #the error vector, for the Kinematic Block
 errors=[error_int,h,h_p]
 #the kinematic block, (Simple PID type controller for 2nd order system) produce the Desired Acceleration.
 qDes_pp[i]=KB(3,Λ,errors,qN_pp[i])
 # If we turn RFPT on, and from the 10th cycle (make sure, the buffers is full)
 if Adaptive==1 && i>10
   #Deformate the system
   qDef_pp[i]=G_SISO(past_input,past_response,qDes_pp[i],K,B,A)
 else
   #if RFPT is off, pass the desired value (thats why it is a simple PID controller)
   qDef_pp[i]=qDes_pp[i]
 end #if
 #store the output of G function, for the next cycle
 past_input=qDef_pp[i]
 # create the control signal
 u[i]=approx(q[i],q_p[i],qDef_pp[i])
 # get the realized response from our system
 q_pp[i]=exact(q[i],q_p[i],u[i])
 # store it for the next cycle (G function needs it.)
 past_response=q_pp[i]
 # Integral of the errror  for the next cycle
 error_int=Integ(error_int,δt,h)
 # Integrals to produce realized q_p and q
 q_p[i+1]=Integ(q_p[i],δt,q_pp[i])
 q[i+1]=Integ(q[i],δt,q_p[i])
end #for
#the last step
t[LONG]=δt*LONG
qN[LONG]=Ampl*sin(ω*t[LONG])
qN_p[LONG]=Ampl*ω*cos(ω*t[LONG])
qN_pp[LONG]=-Ampl*ω^2*sin(ω*t[LONG])
h=qN[LONG]-q[LONG]
h_p=qN_p[LONG]-q_p[LONG]
errors=[error_int,h,h_p]
qDes_pp[LONG]=KB(3,Λ,errors,qN_pp[LONG])
if Adaptive==1
  qDef_pp[LONG]=G_SISO(past_input,past_response,qDes_pp[LONG],K,B,A)
else
  qDef_pp[LONG]=qDes_pp[LONG]
end #if
u[LONG]=approx(q[LONG],q_p[LONG],qDef_pp[LONG])
q_pp[LONG]=exact(q[LONG],q_p[LONG],u[LONG])

#figures

using PyPlot
figure(1)
grid("on")
title("Nominal and realized traj. Nominal: red Realized: blue")
plot(t,qN,color="red")
plot(t,q,color="blue","r--")

figure(2)
grid("on")
title("the control signal")
plot(t,u,color="orange")

figure(3)
grid("on")
title("Tracking Error")
plot(t,qN-q,color="blue")

figure(4)
grid("on")
title("Trajectory Phase Space")
plot(qN,qN_p,color="red")
plot(q,q_p,color="blue","r--")

figure(5)
grid("on")
title("Accelerations: Nominal: red, Desired: orange, Realized: Blue")
plot(t,qN_pp,color="red")
plot(t,qDes_pp,color="orange","r-.")
plot(t,q_pp,color="blue","r--")
