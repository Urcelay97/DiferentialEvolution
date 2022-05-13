include("testfunctions.jl")
using DiferentialEvolution
using Test

@testset "DiferentialEvolution.jl" begin

    @testset "Classic Diferential Evolution" begin

        @testset "2-Dimensioanl tests" begin
            AckleyT = DE_classic(Ackley,1000,300, [-32.768,-32.768], [32.768,32.768],0.9,0.5)
            SchwefelT = DE_classic(Schwefel,1000,300, [-500,-500], [500,500],0.9,0.5)
            GriewankT = DE_classic(Griewank,1000,300, [-600,-600], [600,600],0.9,0.5)
            RosenbrockT = DE_classic(Rosenbrock,1000,300, [-5,-5], [10,10],0.9,0.5)
            MichalewiczT = DE_classic(Michalewicz,1000,300, [0,0], [pi,pi],0.9,0.5)

            #Ackley     critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(AckleyT[1],0, atol=1e-8)
            #Schwefel   critic value = -418.9829 d   ///  parameters = [420.9687,420.9687,420.9687,420.9687,420.9687]
            D = 2
            @test isapprox(SchwefelT[1],-418.9829*D, atol=1e-4)
            #Griewank   critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(GriewankT[1],0, atol=1e-8)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(RosenbrockT[1],0, atol=1e-8)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(MichalewiczT[1],-1.8013, atol=1e-4)
        end


        @testset "5-Dimensional tests" begin
            D = 5
            AckleyT = DE_classic(Ackley,1000,300, [-32.768 for i in 1:D], [32.768 for i in 1:D],0.9,0.5)
            SchwefelT = DE_classic(Schwefel,1000,300, [-500 for i in 1:D], [500 for i in 1:D],0.9,0.5)
            GriewankT = DE_classic(Griewank,1000,300, [-600 for i in 1:D], [600 for i in 1:D],0.9,0.5)
            RosenbrockT = DE_classic(Rosenbrock,1500,300, [-5 for i in 1:D], [10 for i in 1:D],0.9,0.5)
            MichalewiczT = DE_classic(Michalewicz,1000,300, [0 for i in 1:D], [pi for i in 1:D],0.9,0.5)

            #Ackley     critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(AckleyT[1],0, atol=1e-8)
            #Schwefel   critic value = -418.9829 d   ///  parameters = [420.9687,420.9687,420.9687,420.9687,420.9687]
            @test isapprox(SchwefelT[1],-418.9829*D, atol=1e-4)
            #Griewank   critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(GriewankT[1],0, atol=1e-8)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(RosenbrockT[1],0, atol=1e-8)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(MichalewiczT[1],-4.687658, atol=1e-6)
        end

        @testset "15 - Dimensional tests" begin
            D = 15
            AckleyT = DE_classic(Ackley,3000,700, [-32.768 for i in 1:D], [32.768 for i in 1:D],0.9,0.5)
            SchwefelT = DE_classic(Schwefel,3000,700, [-500 for i in 1:D], [500 for i in 1:D],0.9,0.5)
            GriewankT = DE_dither(Griewank,3000,700, [-600 for i in 1:D], [600 for i in 1:D],0.9,0.5)
            RosenbrockT = DE_dither(Rosenbrock,3000,700, [-5 for i in 1:D], [10 for i in 1:D],0.9,0.5)
            #MichalewiczT = DE_classic(Michalewicz,1000,300, [0 for i in 1:D], [pi for i in 1:D],0.9,0.5)

            #Ackley     critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(AckleyT[1],0, atol=1e-2)
            #Schwefel   critic value = -418.9829 d   ///  parameters = [420.9687,420.9687,420.9687,420.9687,420.9687]
            @test isapprox(SchwefelT[1],-418.9829*D, atol=1e-2)
            #Griewank   critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(GriewankT[1],0, atol=1e-2)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(RosenbrockT[1],0, atol=1e-2)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            #@test isapprox(MichalewiczT[1],-12.952351, atol=1e-6)
        end
#=
        @testset "30 - Dimensional tests" begin
            D = 30
            AckleyT = DE_classic(Ackley,5000,700, [-32.768 for i in 1:D], [32.768 for i in 1:D],0.9,0.5)
            SchwefelT = DE_classic(Schwefel,5000,700, [-500 for i in 1:D], [500 for i in 1:D],0.9,0.5)
            GriewankT = DE_classic(Griewank,5000,700, [-600 for i in 1:D], [600 for i in 1:D],0.9,0.5)
            RosenbrockT = DE_classic(Rosenbrock,5000,700, [-5 for i in 1:D], [10 for i in 1:D],0.9,0.5)
            #MichalewiczT = DE_classic(Michalewicz,1000,300, [0 for i in 1:D], [pi for i in 1:D],0.9,0.5)

            #Ackley     critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(AckleyT[1],0, atol=1e-2)
            #Schwefel   critic value = -418.9829 d   ///  parameters = [420.9687,420.9687,420.9687,420.9687,420.9687]
            @test isapprox(SchwefelT[1],-418.9829*D, atol=1e-2)
            #Griewank   critic value = 0   ///  parameters = [0,0,0,0,0]
            @test isapprox(GriewankT[1],0, atol=1e-2)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            @test isapprox(RosenbrockT[1],0, atol=1e-2)
            #Rosenbrock critic value = 0   ///  parameters = [1,1,1,1,1]
            #@test isapprox(MichalewiczT[1],-4.687658, atol=1e-6)
        end
    =#  
    end

    @testset "Dither Diferential Evolution" begin
    
    end

    @testset "Jither Diferential Evolution" begin
        
    end
    
end
