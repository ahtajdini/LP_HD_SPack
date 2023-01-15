
module sphere_main

include("./seed_initializer.jl")
include("./loss.jl")
include("./guess_functions.jl")

using .Loss, .Guess, .rSeed


hyper = Dict("precision"=>250,"threshold" => 10^(-25),
"order" => 5, "lr" => 1, "nmax_ph1" =>50, "nmax_ph2"=>100, "nmax_ph3"=> 200,
"step_ph0"=>1,"step_ph1" => 5, "step_ph2" => 10, "step_ph3" =>30,
 "bet_gues_diff" => 25 )


setprecision(BigFloat, hyper["precision"],base=10)


function Block1(cent,nmax,seed,res_dic)

    n= length(seed)
    pre_ans = [seed;-ones(BigFloat,n)]
    res_dic[n] = Loss.newtonoptimizer(cent, pre_ans, hyper["lr"], hyper["threshold"] );

    for i in range(n+1,nmax)
        pre_root = Guess.smartGuess(res_dic[i-1][1:i-1,1],i)
        pre_coeff=  -ones(BigFloat,i)
        pre_ans= [ pre_root ; pre_coeff]
        res_dic[i] = Loss.newtonoptimizer(cent,pre_ans,hyper["lr"], hyper["threshold"])

    end
end

function Block2(cent,nmin,nmax,step, res_dic,tab)

    for i in range(nmin,nmax,step=step)

        pre_root = Guess.better_guess2(tab, res_dic, i )
        pre_coeff=  -ones(BigFloat,i)
        pre_ans= [ pre_root ; pre_coeff]

        res_dic[i] = Loss.newtonoptimizer(cent,pre_ans, hyper["lr"], hyper["threshold"])

        push!(tab,i)
    end
end

#cent is central charge==D/2, 
#N is the number of roots. seed is the initial values for roots,
#  it is often either the first three or first four roots
function lpSpacking(D,N,res_dic)
    cent =BigFloat(D)/2
    
    seed = rSeed.seed(D)

    Block1(cent,hyper["nmax_ph1"],seed,res_dic)


    if N>hyper["nmax_ph1"] 

        tab = collect(hyper["nmax_ph1"]-hyper["bet_gues_diff"]:hyper["nmax_ph1"])
        Block2(cent, hyper["nmax_ph1"]+hyper["step_ph1"],min(N,hyper["nmax_ph2"]),hyper["step_ph1"], res_dic,tab)
    end
    if N>hyper["nmax_ph2"]
        Block2(cent, hyper["nmax_ph2"]+hyper["step_ph2"],min(N,hyper["nmax_ph3"]),hyper["step_ph2"],res_dic,tab)
    end

    if N> hyper["nmax_ph3"]

        Block2(cent,hyper["nmax_ph3"]+hyper["step_ph3"], N,hyper["step_ph3"] , res_dic,tab)
    end

end

end






