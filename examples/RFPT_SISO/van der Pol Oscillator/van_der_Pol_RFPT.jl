#############################################
# Van der Pol Oscillator Controlled by RFPT #
#############################################
Adaptive=1
######################
# Control Parameters #
######################
K=1e5
B=-1
A=1e-5

#########
# Time  #
#########
δt=1e-3
LONG=Int(2e4)

########################################
# Kinematic Block Parameter (2nd order)#
########################################
Λ=3

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

############################
# Functions for simulation #
############################

# function nom_traj(t)
#   #egys=[1,-1.5,2]
#   qN=Ampl*sin(ω*δt*t)
#   qN_p=Ampl*ω*cos(ω*δt*t)
#   qN_pp=-Ampl*ω^2*sin(ω*δt*t)
# 	#println(qN[t,:])
#   return [qN;qN_p;qN_pp]
# end

#######################
# The Kinematic Block #
#######################

# function Kinematic_Block(t,err_int)
# 	#Compute the Error
# 	err=(qN-q)#'
# 	#Compute the 1st time derivative of the Error
# 	err_p=(qN_p-q_p)#'
# 	#Compute te integral of the Error (with Euler formula)
#    #println(err)
# 	err_i=err_int+err*δt
#   #println(q,q_p,q_pp)
# 	#Compute the Desired acceleration (2nd time derivative)
# 	#out=qN_pp[t]'+3*Λ^2*err+3*Λ*err_p+Λ^3*err_i
#   out=qN_pp+3*Λ^2*err+3*Λ*err_p+Λ^3*err_i #Has to be colum vector
#   #println(out)
# 	#println(err_i)
#    #println(typeof(err_i))
# 	return [out err_i]
# end

###############
# Exact Model #
###############
function Exact(u,q,q_p)
  q_pp=(u+μₑ*(1-q^2)*q_p-ωₑ^2*q-αₑ*q^3-λₑ^5)/mₑ
  return q_pp
end

################
# Approx Model #
################

function Approx(q,q_p,q_pp)
  u=mₐ*q_pp-μₐ*(1-q^2)*q_p+ωₐ^2*q+αₐ*q^3+λₐ^5
  return u
end

function sigmoid(x)
	s=x/(1+abs(x))
  return s
end

function G(past_input,past_response,xDnow)#(f,x,xᵈ)
  out=(K+past_input) * (1 + B * sigmoid(A * (past_response - xDnow))) - K
  return out #q_ppReq
end

##############################
# Define arrays for Plotting #
##############################
#time
time_mem=zeros(LONG)
#nominal trajectory
#qN_mem=zeros(LONG)
#qN_p_mem=zeros(LONG)
#qN_pp_mem=zeros(LONG)

#q_ppDes_mem=zeros(LONG)
#x=0
#f=0
#q_ppDef_mem=zeros(LONG)
#Q_mem=zeros(LONG)
#err_int=0
#qN=0
#qN_p=0
#qN_pp=0
#q=0
#q_p=0
#q_pp=0
#q_mem=zeros(LONG)
##################
# The simulation #
##################
q_mem=zeros(LONG) #A megvalósult pálya
q_p_mem=zeros(LONG)

qN_mem=zeros(LONG) #A nominális pálya
qN_p_mem=zeros(LONG)

u_mem=zeros(LONG) #A szanályozó jel elmentése

qN_pp_mem=zeros(LONG) #A nominális gyorsulások értéke
q_pp_Des_mem=zeros(LONG)#A PID korrekciós adatok
q_pp_Req_mem=zeros(LONG) #Az Adaptívan torzított jelre
q_pp_mem=zeros(LONG) #Az Adaptívan torzított jelre

h_int=0 #A követési hiba integrálja

past_input=0
past_response=0

for t=1:LONG-1
  #store time for plotting
  time_mem[t]=t*δt
  #define the nominal trajectory
  qN_mem[t]=Ampl*sin(ω*δt*t)
  qN_p_mem[t]=ω*Ampl*cos(ω*δt*t)
  qN_pp=-ω^2*Ampl*sin(ω*δt*t)
  qN_pp_mem[t]=qN_pp
  #the kinematic
  q_ppDes=Λ^3*h_int+3*Λ^2*(qN_mem[t]-q_mem[t])+3*Λ*(qN_p_mem[t]-q_p_mem[t])+qN_pp
  q_pp_Des_mem[t]=q_ppDes
  if Adaptive==1 && t>10
    q_ppReq=G(past_input,past_response,q_ppDes)
  else
    q_ppReq=q_ppDes
  end
  q_pp_Req_mem[t]=q_ppReq
  u=Approx(q_mem[t],q_p_mem[t],q_ppReq)
  u_mem[t]=u
  #println(Q)
  past_input=q_ppReq

  q_pp=Exact(u,q_mem[t],q_p_mem[t])
  #println(x)
  #println(x)
  q_pp_mem[t]=q_pp

  past_response=q_pp
  q_p_mem[t+1]=q_p_mem[t]+δt*q_pp
  q_mem[t+1]=q_mem[t]+δt*q_p_mem[t]
  h_int=h_int+δt*(qN_mem[t]-q_mem[t])

end #for t=1:LONG-1

using(PyPlot)
figure(1)
grid("on")
title("Tracking Error vs Time")
xlabel("Time [s]")
ylabel("Tracking Error")
plot(time_mem[11:LONG-1],qN_mem[11:LONG-1]-q_mem[11:LONG-1])

figure(2)
grid("on")
title("Nominal and Realized Trajectory")
xlabel("Time [s]")
ylabel("Nominal and Realized Trajectory")
plot(time_mem[3:LONG-1],qN_mem[3:LONG-1],color="red")
plot(time_mem[3:LONG-1],q_mem[3:LONG-1],color="green")
