include("testfunctions.jl")
using DiferentialEvolution
using PlotlyJS

#Plots
len=1000
x = LinRange(0,pi,len)
y = LinRange(0,pi,len)
X = ones(len)*x'
Y = y*ones(len)'
Z = Michalewicz.(X,Y)

critic = DE_jither(Michalewicz,1000,500, [0,0], [pi,pi],0.9,0.5)
p_init = contour(x =x, y= y,z = Z)
crit = scatter(x=[critic[2][1]],y=[critic[2][2]], mode="markers",color=2)
plot([p_init,crit])
