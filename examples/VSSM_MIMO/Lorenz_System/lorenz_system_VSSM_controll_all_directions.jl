############################
# Lorenz System  (MIMO)    #
# by VSSM (all directions) #
############################

###########################################
# Controll parameters for VSSM controller #
###########################################

K=500
w=1

#############################################
# Parameters to design a Nominal Trajectory #
#############################################

A₁=2
ω₁=0.5
A₂=3.8
ω₂=0.7
A₃=1.5
ω₃=0.9

#################
# Time variable #
#################

# cycle time in seconds
δt=1e-3
# The lenght of the simulation in seconds
LONG=Int(2e4)

Λ=12

##########################
# Exact model parameters #
##########################
σₑ=5
βₑ=8/3
ρₑ=40

############################
# Approx. model parameters #
############################
σₐ=5
βₐ=8/3
ρₐ=40

#Arrays

xN=zeros(LONG)
xN_p=zeros(LONG)
yN=zeros(LONG)
yN_p=zeros(LONG)
zN=zeros(LONG)
zN_p=zeros(LONG)

x_p_des=zeros(LONG)
y_p_des=zeros(LONG)
z_p_des=zeros(LONG)

u_x=zeros(LONG)
u_y=zeros(LONG)
u_z=zeros(LONG)

t=zeros(LONG)

S_x=zeros(LONG)
S_p_x=zeros(LONG)
S_y=zeros(LONG)
S_p_y=zeros(LONG)
S_z=zeros(LONG)
S_p_z=zeros(LONG)

x=zeros(LONG)
y=zeros(LONG)
z=zeros(LONG)
x_p=zeros(LONG)
y_p=zeros(LONG)
z_p=zeros(LONG)

#initial conditions
hint_x=0
hint_y=0
hint_z=0

x[1]=A₁*sin(ω₁*δt)
y[1]=A₂*sin(ω₂*δt)
z[1]=A₃*sin(ω₃*δt)


##############
# SIMULATION #
##############

for i=1:LONG-1 #i counts the cycles of the simulation
  t[i]=δt*i
  #create Nominal Trajectory
  xN[i]=A₁*sin(ω₁*t[i])
  xN_p[i]=A₁*ω₁*cos(ω₁*t[i])
  #xN_p[i]=-A₁*ω₁^2*sin(ω₁*t[i])

  yN[i]=A₂*sin(ω₂*t[i])
  yN_p[i]=A₂*ω₂*cos(ω₂*t[i])
  #yN_p[i]=-A₂*ω₂^2*sin(ω₂*t[i])

  zN[i]=A₃*sin(ω₃*t[i])
  zN_p[i]=A₃*ω₃*cos(ω₃*t[i])
  #zN_p[i]=-A₃*ω₃^2*sin(ω₃*t[i])
  # Compute  error (h) and error integral (hint) and the first derivative of the error (h_p)
  h_x=xN[i]-x[i]
  h_y=yN[i]-y[i]
  h_z=zN[i]-z[i]

  hint_x=hint_x+δt*h_x
  hint_y=hint_y+δt*h_y
  hint_z=hint_z+δt*h_z

  h_p_x=xN_p[i]-x_p[i]
  h_p_y=yN_p[i]-y_p[i]
  h_p_z=zN_p[i]-z_p[i]

  #Error metric (S), and the first derivative of the Error metric (S_p)
  S_x[i]=Λ*hint_x+h_x
  S_p_x[i]=Λ*h_x+h_p_x
  S_y[i]=Λ*hint_y+h_y
  S_p_y[i]=Λ*h_y+h_p_y
  S_z[i]=Λ*hint_z+h_z
  S_p_z[i]=Λ*h_z+h_p_z
  #Kinematic Block (PID type feedback)
  x_p_des[i]=xN_p[i]+Λ*h_x+K*tanh(S_x[i]/w)
  y_p_des[i]=yN_p[i]+Λ*h_y+K*tanh(S_y[i]/w)
  z_p_des[i]=zN_p[i]+Λ*h_z+K*tanh(S_z[i]/w)
  # Create the control signal
  u_x[i]=x_p_des[i]-σₐ*(y[i]-x[i])
  u_y[i]=y_p_des[i]-x[i]*(ρₐ-z[i])
  u_z[i]=z_p_des[i]-y[i]+βₐ*z[i]
  # Get the system "reaction" (the exact model of the system)
  x_p[i]=σₑ*(y[i]-x[i])+u_x[i]
  y_p[i]=x[i]*(ρₑ-z[i])-y[i]+u_y[i]
  z_p[i]=x[i]*y[i]-βₑ*z[i]+u_z[i]
  # Integrate back to x, y, z
  x[i+1]=x[i]+δt*x_p[i]
  y[i+1]=y[i]+δt*y_p[i]
  z[i+1]=z[i]+δt*z_p[i]
end #for

using PyPlot

figure(1)
grid("on")
title("Nominal and Realized trajectorys")
xlabel(L"Time $[s]$")
ylabel(L"q $[m]$")
plot(t[3:LONG-1],xN[3:LONG-1],color="#7684FF")
plot(t[3:LONG-1],yN[3:LONG-1],color="#2A40FF")
plot(t[3:LONG-1],zN[3:LONG-1],color="#3B427F")

plot(t[3:LONG-1],x[3:LONG-1],color="#FFAA41","r--")
plot(t[3:LONG-1],y[3:LONG-1],color="#FF9F28","r--")
plot(t[3:LONG-1],z[3:LONG-1],color="#B2690E","r--")

figure(2)
grid("on")
title("Tracking Error")
xlabel(L"Time $[s]$")
ylabel(L"error $[m]$")
plot(t[3:LONG-1],xN[3:LONG-1]-x[3:LONG-1],color="red")
plot(t[3:LONG-1],yN[3:LONG-1]-y[3:LONG-1],color="green")
plot(t[3:LONG-1],zN[3:LONG-1]-z[3:LONG-1],color="blue")

figure(3)
grid("on")
title("Phase Space of Error Metric (S_p vs S)")
plot(S_x[3:LONG-1],S_p_x[3:LONG-1],color="red")
plot(S_y[3:LONG-1],S_p_y[3:LONG-1],color="green","r--")
plot(S_z[3:LONG-1],S_p_z[3:LONG-1],color="blue","r-.")

figure(4)
grid("on")
title("Phase Space of motion (x_p vs x)")
plot(xN[1:LONG-1],xN_p[1:LONG-1],color="red")
plot(x[1:LONG-1],x_p[1:LONG-1],color="blue","r--")
plot(yN[1:LONG-1],yN_p[1:LONG-1],color="#7684FF")
plot(y[1:LONG-1],y_p[1:LONG-1],color="#FFAA41","r--")
plot(zN[1:LONG-1],zN_p[1:LONG-1],color="#2A40FF")
plot(z[1:LONG-1],z_p[1:LONG-1],color="#FF9F28","r-.")
