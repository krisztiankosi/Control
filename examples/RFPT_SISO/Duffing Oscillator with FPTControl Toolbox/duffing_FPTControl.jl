
 #Duffing Oscillator
 #Controller: RFPT (SISO)
 #using FPTControl toolbox

using FPTControl

Adaptive=1
LONG=Int(1e4)
δt=1e-3

#Kinematic Parameter
Λ=3

#Nominal Trajectory
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

# Approx model parameters
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

################
# Approx Model #
################

function approx(q,q_p,q_pp)
  out=q_pp#-αₐ*q-δₐ*q_p-βₐ*q^3
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

error_int=0
past_input=0
past_response=0

for i=1:LONG-1
  t[i]=i*δt
  #Nominal trajectory
  qN[i]=Ampl*sin(ω*t[i])#sin(ω*time_mem[t])
  qN_p[i]=ω*Ampl*cos(ω*t[i])
  qN_pp[i]=-ω^2*Ampl*sin(ω*t[i])
  #Kinematic Block

  qDes_pp[i]=Λ^3*error_int+3*Λ^2*(qN[i]-q[i])+3*Λ*(qN_p[i]-q_p[i])+qN_pp[i]
  qDes_pp[i]
  # Deformation
  if Adaptive==1 && i>10
    qDef_pp[i]=G_SISO(past_input,past_response,qDes_pp[i],K,B,A)
  #  println(qDef_pp[i])
  else
    qDef_pp[i]=qDes_pp[i]
  end
  #Approx model
  u[i]=approx(q[i],q_p[i],qDef_pp[i])
  #u=qDef_pp_mem[t]
  past_input=qDef_pp[i]
  #Exact model
  q_pp[i]=exact(q[i],q_p[i],u[i])
  past_response=q_pp[i]
  #Integrals
  #q_p[i+1]=q_p[i]+δt*q_pp[i]
   q_p[i+1]=Integ(q_p[i],δt,q_pp[i])
  #q[i+1]=q[i]+δt*q_p[i]
  q[i+1]=Integ(q[i],δt,q_p[i])
  error_int=Integ(error_int,δt,qN[i]-q[i])
  #error_int=error_int+δt*(qN[i]-q[i])
end #for (simulation)

using PyPlot

figure(1)
title("Tracking error")
plot(t[3:LONG-1],qN[3:LONG-1]-q[3:LONG-1])
grid("on")
xlabel("Time [s]")
ylabel("error [m]")
#xlabel(L"\frac{x}{y}")   ##(L"") latex-es felirathoz

figure(2)
title("Nominal and Simulated")
plot(t[3:LONG-1],qN[3:LONG-1],color="red")
plot(t[3:LONG-1],q[3:LONG-1],color="green","r--")
grid("on")
xlabel("Time [s]")
ylabel("q [m]")
