include("testfunctions.jl")
using DiferentialEvolution
using PlotlyJS

#Plots
len=500
x = LinRange(-600,600,len)
y = LinRange(-600,600,len)
X = ones(len)*x'
Y = y*ones(len)'
Z = Eggholder.(X,Y)
Eggholder(-2.5153506020027042e64, 1.79281509087981e64)
@time critic1 = DE_JADE(Eggholder, 500, 25, [400,400], [600,600], 0.2, 1/20)
@time critic2 = DE_dither(Eggholder, 1000, 500, [-512,-512], [512,512],0.9,0.5)
p_init = contour(x =x, y= y,z = Z)
crit = scatter(x=[critic[2][1]],y=[critic[2][2]], mode="markers",color=2)
plot([p_init,crit])
