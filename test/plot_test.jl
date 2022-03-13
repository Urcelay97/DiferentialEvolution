include("testfunctions.jl")
using DiferentialEvolution
using PlotlyJS


#Plots
len=1000
x = LinRange(-5.12,5.12,len)
y = LinRange(-5.12,5.12,len)
X = ones(len)*x'
Y = y*ones(len)'
Z = Drop_wave.(X,Y)

critic = DE_jither(Drop_wave,5000,500, [-10,-10], [10,10],0.9,0.5)
p_init = contour(x =x, y= y,z = Z)
crit = scatter(x=[critic[2][1]],y=[critic[2][2]], mode="markers",color=2)
plot([p_init,crit])
