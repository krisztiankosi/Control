######################
# Duffing Oscillator #
######################

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
  out=q_pp-αₐ*q-δₐ*q_p-βₐ*q^3
  return out #u
end

########################
# Deformation Function #
########################
function sigmoid(x)
	s=x/(1+abs(x))
  return s
end
function G(past_input,past_response,desired)
  out=(K+past_input)*(1 + B * tanh(A * (past_response - desired))) - K
  return out
end #qDef_pp

######################################
# Arrays for Plotting and simulation #
######################################

qN_mem=zeros(LONG)
qN_p_mem=zeros(LONG)
qN_pp_mem=zeros(LONG)

qDes_pp_mem=zeros(LONG)

qDef_pp_mem=zeros(LONG)
u_mem=zeros(LONG)

q_mem=zeros(LONG)
q_p_mem=zeros(LONG)
q_pp_mem=zeros(LONG)

time_mem=zeros(LONG)

error_int=0
past_input=0
past_response=0

#simulation
for t=1:LONG-1
  time_mem[t]=t*δt
  #Nominal trajectory
  qN_mem[t]=Ampl*sin(ω*δt*t)#sin(ω*time_mem[t])
  qN_p_mem[t]=ω*Ampl*cos(ω*δt*t)
  qN_pp_mem[t]=-ω^2*Ampl*sin(ω*δt*t)
  #Kinematic Block
  qDes_pp_mem[t]=Λ^3*error_int+3*Λ^2*(qN_mem[t]-q_mem[t])+3*Λ*(qN_p_mem[t]-q_p_mem[t])+qN_pp_mem[t]
  # Deformation
  if Adaptive==1 && t>10
    qDef_pp_mem[t]=G(past_input,past_response,qDes_pp_mem[t])
  else
    qDef_pp_mem[t]=qDes_pp_mem[t]
  end
  #Approx model
  u=approx(q_mem[t],q_p_mem[t],qDef_pp_mem[t])
  u_mem[t]=u
  past_input=qDef_pp_mem[t]
  #Exact model
  q_pp_mem[t]=exact(q_mem[t],q_p_mem[t],u)
  past_response=q_pp_mem[t]
  #Integrals
  q_p_mem[t+1]=q_p_mem[t]+δt*q_pp_mem[t]
  q_mem[t+1]=q_mem[t]+δt*q_p_mem[t]
  error_int=error_int+δt*(qN_mem[t]-q_mem[t])
end #for (simulation)

using PyPlot

figure(1)
title("Tracking error")
plot(time_mem[3:LONG-1],qN_mem[3:LONG-1]-q_mem[3:LONG-1])
grid("on")
xlabel("Time [s]")
ylabel("error [m]")
#xlabel(L"\frac{x}{y}")   ##(L"") latex-es felirathoz

figure(2)
title("Nominal and Simulated")
plot(time_mem[3:LONG-1],qN_mem[3:LONG-1],color="red")
plot(time_mem[3:LONG-1],q_mem[3:LONG-1],color="green")
grid("on")
xlabel("Time [s]")
ylabel("q [m]")
