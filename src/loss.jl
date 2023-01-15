

module Loss


using LinearAlgebra, Dates, Printf



function vac_vector(vac,n,N1)

    aux = ones(BigFloat,(4*n,1))
    aux[2,1] = N1-2 * vac
    
    for i in range(3,4*n)
        aux[i,1] = ((-i+3-N1).* aux[i-2,1]+(2*i+N1-4-2*vac) .* aux[i-1,1])./ (i-1)
    end
    return aux
end

function Laguerre_root(root_list,N1)
    n = length(root_list)
    row = ones(BigFloat,(4*n,n))
    row[2,:] = N1 .- 2 .* root_list
    for i in 3:4*n
        row[i,:] = (-(N1-3+i) .* row[i-2,:] + (2*i+N1-4 .- 2 .* root_list) .* row[i-1,:] ) ./ (i-1)
    end
    return row
end

function derivative_Laguerre(rl,N1,rows)
    n = length(rl)
    drow = zeros(BigFloat,(4*n,n))
    drow[2,:] .= -2
    for i in 3:4*n
        drow[i,:] = ( (i-2 .- 2 * rl) .* drow[i-1,:] - 2*(i + N1-3) .* rows[i-2,:]) ./ (i-2)
    end
    return drow
end

#Goal: set loss=0. vac is the value
function obj(vac,full_list,N1)

    n=Int(length(full_list)/2)

    root_list= full_list[1:n]
    coef_list = full_list[n+1:2*n]

    one_vec = ones(BigFloat,(2*n,1))

    aux= vac_vector(vac,n,N1)

    row = Laguerre_root(root_list,N1)

    column= transpose(row[2:2:4*n,:] ./ aux[2:2:4*n,:])
    
    maxim = column[findmax(abs.(column),dims=2)[2]]
    mat = column./maxim
    result =transpose(mat)*coef_list+ one_vec
    return result
end


function fT(cent,X)
    return obj(0, X, cent)
    
end


function Jacobian(vac,full_list,N1)
    n=Int(length(full_list)/2)
    root_list= full_list[1:n]
    coef_list = full_list[n+1:2*n]
    aux = vac_vector(vac,n,N1)
    row = Laguerre_root(root_list,N1)
    drow = derivative_Laguerre(root_list,N1,row)
    seed = transpose(row[2:2:4*n,:] ./ aux[2:2:4*n,:])
    dseed = transpose(drow[2:2:4*n,:] ./ aux[2:2:4*n,:])
    _, pos = findmax(abs.(seed),dims=2)
    maxim = seed[pos]
    mat = seed ./ maxim
    dmat = dseed ./ maxim - (seed .* dseed[pos]) ./ (maxim .* maxim)
    return transpose(vcat(dmat .* coef_list, mat))
end

function Df(cent,X) 
    return Jacobian(0,X,cent)
end


function newtonoptimizer(cent,X0; thresh=10^(-25))
    MAX_NUM= 10^(20)
    t0=now()
    n= Int(length(X0)/2)
    X=copy(X0)
    f= fT(cent,X)

    cond = true
    count=0
    loss = (f'f)[1]
    stopreason ="Undetermined"
    while cond
        jac = Df(cent,X)
        count += 1
        dX = (-inv(jac))*f
        mu=1
        X += mu*dX
        f = fT(cent,X)
        loss = (f'f)[1]
        if loss < thresh
            cond = false
            stopreason = "Success"
        elseif count>100
            stopreason = "Max Iteration"
            cond=false
        elseif loss >  MAX_NUM
            stopreason ="MAX_NUM"
            cond = false
        end
    end

    ntime =  Dates.canonicalize(Dates.CompoundPeriod(now()-t0))
    @printf("%s | N = %d | error = %.2e | iters = %d | time =  %s \n", stopreason,n,loss,count, ntime)
    return X
end


end
