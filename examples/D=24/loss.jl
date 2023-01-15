
"""This file includes the objective function, jacobian and optimizer which 
implements newton method"""


module Loss

using LinearAlgebra, Dates, Printf


#this function computes the Laguerre polynomials evaluated at origin
function vacuum_vector(vac,n,N1)

    aux = ones(BigFloat,(4*n,1))
    aux[2,1] = N1-2 * vac
    
    for i in range(3,4*n)
        aux[i,1] = ((-i+3-N1).* aux[i-2,1]+(2*i+N1-4-2*vac) .* aux[i-1,1])./ (i-1)
    end
    return aux
end


#This computes the laguerre polynomials for the candidate list of roots
function Laguerre_list(root_list,N1)
    n = length(root_list)
    row = ones(BigFloat,(4*n,n))
    row[2,:] = N1 .- 2 .* root_list
    for i in 3:4*n
        row[i,:] = (-(N1-3+i) .* row[i-2,:] + (2*i+N1-4 .- 2 .* root_list) .* row[i-1,:] ) ./ (i-1)
    end
    return row
end


#This computes the derivate of Laguerre polys using recurcive relations.
function derivative_Laguerre(rl,N1,rows)
    n = length(rl)
    drow = zeros(BigFloat,(4*n,n))
    drow[2,:] .= -2
    for i in 3:4*n
        drow[i,:] = ( (i-2 .- 2 * rl) .* drow[i-1,:] - 2*(i + N1-3) .* rows[i-2,:]) ./ (i-2)
    end
    return drow
end

#The objective function. Inputs are value of vacuum (always 0 for sphere packign)
#list of roots and coefficients, and parameter N1 which is equal to cent for sphere packing
function obj(vac,full_list,N1)

    n=Int(length(full_list)/2)

    root_list= full_list[1:n]
    coef_list = full_list[n+1:2*n]

    one_vec = ones(BigFloat,(2*n,1))

    aux= vacuum_vector(vac,n,N1)

    row = Laguerre_list(root_list,N1)

    column= transpose(row[2:2:4*n,:] ./ aux[2:2:4*n,:])
    
    maxim = column[findmax(abs.(column),dims=2)[2]]
    mat = column./maxim
    result =transpose(mat)*coef_list+ one_vec
    return result
end

#The jacobian of the objective function
function Jacobian(vac,full_list,N1)
    n=Int(length(full_list)/2)
    root_list= full_list[1:n]
    coef_list = full_list[n+1:2*n]
    aux = vacuum_vector(vac,n,N1)
    row = Laguerre_list(root_list,N1)
    drow = derivative_Laguerre(root_list,N1,row)
    seed = transpose(row[2:2:4*n,:] ./ aux[2:2:4*n,:])
    dseed = transpose(drow[2:2:4*n,:] ./ aux[2:2:4*n,:])
    _, pos = findmax(abs.(seed),dims=2)
    maxim = seed[pos]
    mat = seed ./ maxim
    dmat = dseed ./ maxim - (seed .* dseed[pos]) ./ (maxim .* maxim)
    return transpose(vcat(dmat .* coef_list, mat))
end

#Same as objective function. obj function is applicable to broader set
#of problems than sphere packing
function fT(cent,X)
    return obj(0, X, cent)
    
end

#Same as Jacobian
function Df(cent,X) 
    return Jacobian(0,X,cent)
end


#newtonoptimizer function. Inputs are cent defined as D/2,
# X0 list of roots (Deltas) and coeffs (d_Delta), lr learning rate, thresh is threshold in loss
function newtonoptimizer(cent,X0, lr, thresh)
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

        X += lr*dX
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