include("testfunctions.jl")
using DiferentialEvolution
using PlotlyJS

#Plots
#=
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
=#

f(x) = sinc(x)
x = -5:0.1:5
y = f.(x)
g(par...) = par*[@. x^i for i in 1:length(par)]
error(a,b,c,d,e,f,h,i,k,l,m,n) = sum(abs.(y-g(a,b,c,d,e,f,h,i,k,l,m,n)))
res = DiferentialEvolution.DE_JADE(error,1000,500,[-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100],[100,100,100,100,100,100,100,100,100,100,100,100],0.15,0.2)

a1 = scatter(x=x, y=y)
a2 = scatter(x=x, y=g(res[2]...))
plot([a1,a2])