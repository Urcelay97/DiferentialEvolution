include("testfunctions.jl")
using DiferentialEvolution
using Test
MichalewiczT = DE_classic(Michalewicz,1000,300, [0,0,0,0,0], [pi,pi,pi,pi,pi],0.9,0.5)
@testset "DiferentialEvolution.jl" begin

    @testset "Classic Diferential Evolution" begin

        @testset "5-Dimensional tests" begin

            AckleyT = DE_classic(Ackley,1000,300, [-32.768,-32.768,-32.768,-32.768,-32.768], [32.768,32.768,32.768,32.768,32.768],0.9,0.5)
            SchwefelT = DE_classic(Schwefel,1000,300, [-500,-500,-500,-500,-500], [500,500,500,500,500],0.9,0.5)
            GriewankT = DE_classic(Griewank,1000,300, [-600,-600,-600,-600,-600], [600,600,600,600,600],0.9,0.5)
            RosenbrockT = DE_classic(Rosenbrock,1000,300, [-5,-5,-5,-5,-5], [10,10,10,10,10],0.9,0.5)
            MichalewiczT = DE_classic(Michalewicz,1000,300, [0,0,0,0,0], [pi,pi,pi,pi,pi],0.9,0.5)

            #Ackley     critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(AckleyT[1],0, atol=1e-8)
            #Schwefel   critic value = -418.9829 d   ///  parameters = [420.9687,420.9687,420.9687,420.9687,420.9687]
            d = 5
            @test isapprox(SchwefelT[1],-418.9829*d, atol=1e-4)
            #Griewank   critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(GriewankT[1],0, atol=1e-8)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(GriewankT[1],0, atol=1e-8)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(MichalewiczT[1],-4.687658, atol=1e-6)
        end
        # 2-Dimensional tests

        @testset "2-Dimensioanl tests" begin
        end
    end

    @testset "Dither Diferential Evolution" begin
    
    end

    @testset "Jither Diferential Evolution" begin
        
    end
    
end
