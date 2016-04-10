#############################################
# Van der Pol Oscillator controlled by VSSM #
#############################################

K=500
w=1


#########
# Time  #
#########
δt=1e-3
LONG=Int(2e4)

########################################
# Kinematic Block Parameter (2nd order)#
########################################
Λ=12

####################################
# Parameters for Nominal Trajectory#
####################################

ω=0.5
Ampl=2

##########################
# Exact Model Parameters #
##########################

μₑ=0.4
ωₑ=0.46
αₑ=1
λₑ=0.1
mₑ=1

################################
# Approximate Model Parameters #
################################

μₐ=0.5
ωₐ=0.42
αₐ=0.9
λₐ=0.09
mₐ=0.8

###############
# Exact Model #
###############
function Exact(q,q_p,u)
  q_pp=(u+μₐ*(1-q^2)*q_p-ωₐ^2*q-αₐ*q^3-λₐ^5)/mₑ
  return q_pp
end

################
# Approx Model #
################

function Approx(q,q_p,q_pp)
  u=mₐ*q_pp-μₐ*(1-q^2)*q_p+ωₐ^2*q+αₐ*q^3+λₐ^5
  return u
end



#ErrorMetrics(hint,h,h_p) = Λ^2*hint+2*Λ*h+h_p # "S" block for 2nd order

#KinBlock(S,h,h_p,qN_pp) = K*tanh(S/w)+Λ^2*h+2*Λ*h_p+qN_pp #q_ppDes -t adja vissza

ErrorMetrics(hint,h,h_p) = Λ^2*hint+2*Λ*h+h_p # "S" block for 2nd order

KinBlock(S,h,h_p,qN_pp) = K*tanh(S/w)+Λ^2*h+2*Λ*h_p+qN_pp #q_ppDes -t adja vissza

function nominalTraj(t)
  qN=Ampl*sin(ω*t*δt)
  q_pN=ω*Ampl*cos(ω*t*δt)
  q_ppN=-Ampl*ω^2*sin(ω*t*δt)
  return [qN,q_pN,q_ppN]
end

q_mem=zeros(LONG)
q_p_mem=zeros(LONG)
q_pp_mem=zeros(LONG)
q_ppDes_mem=zeros(LONG)
q_N_mem=zeros(LONG)
q_pN_mem=zeros(LONG)
q_ppN_mem=zeros(LONG)
time_mem=zeros(LONG)
S_mem=zeros(LONG)
u_mem=zeros(LONG)

hint=0
##############
# Simulation #
##############

for t=1:LONG-1
  #create time
  time_mem[t]=δt*t
  #Nominal trajectory for the actual time frame
#  q_temp=nominalTraj(time_mem[t])
  q_temp=nominalTraj(t)
  q_N_mem[t] = q_temp[1]
  q_pN_mem[t] = q_temp[2]
  q_ppN_mem[t] = q_temp[3]
  #errors
  h=q_N_mem[t]-q_mem[t]
  #println(h)
  #Error Integral (Euler integral)
  hint=hint+δt*h
  #The time derivative of the error
  h_p=q_pN_mem[t]-q_p_mem[t]
  #the error metric
  S_mem[t]=ErrorMetrics(hint,h,h_p)
  #the kinematic block
  q_ppDes_mem[t]=KinBlock(S_mem[t],h,h_p,q_ppN_mem[t])
  #println(q_ppDes_mem[t])
  #the aprox model (inverse model)
  u_mem[t]=Approx(q_mem[t],q_p_mem[t],q_ppDes_mem[t])
  #the exact model
  q_pp_mem[t]=Exact(q_mem[t],q_p_mem[t],u_mem[t])
  #	Integrate back with Euler integral
  q_p_mem[t+1]=q_p_mem[t]+δt*q_pp_mem[t]
  q_mem[t+1]=q_mem[t]+q_p_mem[t]*δt
end #for

###############################
# Load matplotlib for Figures #
###############################

using PyPlot
#println(time_mem)
figure(1)
plot(time_mem[3:LONG-1],q_N_mem[3:LONG-1],color="#7684FF")
plot(time_mem[3:LONG-1],q_mem[3:LONG-1],color="#FFAA41")
grid("on")
title("Nominal and Realized Trajectory\n Nominal Trajectory components: bluish\n Realized Trajectory components:orangish ")
xlabel("Time [s]")
ylabel("q [m]")
#plot(time_mem[3:LONG-1],qN[3:LONG-1,1],color="#7684FF")
#plot(time_mem[3:LONG-1],qN[3:LONG-1,2],color="#2A40FF")
#plot(time_mem[3:LONG-1],qN[3:LONG-1,3],color="#3B427F")
#plot(time3d[3:LONG-1],qR[3:LONG-1,:],color="green")
#plot(time_mem[3:LONG-1],qR[3:LONG-1,1],color="#FFAA41")
#plot(time_mem[3:LONG-1],qR[3:LONG-1,2],color="#FF9F28")
#plot(time_mem[3:LONG-1],qR[3:LONG-1,3],color="#B2690E")

figure(2)
plot(time_mem[3:LONG-1],q_N_mem[3:LONG-1]-q_mem[3:LONG-1])
