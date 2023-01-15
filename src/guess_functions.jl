"""This module computes a hot starting for non-linear equations. The equations
are very sensitive to initial guess. Knowing the solutions of equations 
for 1,2,...,N, we come up with a good guess for initial value of N+1 equations.
The idea of better guess function is due to Nima Afkhami-Jeddi """
module Guess

using LinearAlgebra, SpecialPolynomials, Interpolations, Dates 

export smartGuess, better_guess


#newN is in practice =len(oldroots)+1
function smartGuess(oldroots, newN)
    n= length(oldroots)
    rRef = linear_interpolation(range(1,n),oldroots,extrapolation_bc=Line())
    rNew = zeros(BigFloat,(newN,))
    for i in range(1,newN)
        rNew[i] = ((n-newN)*oldroots[1])/n+ newN*rRef(1 + BigFloat((i-1)*(n-1))/BigFloat((newN-1)))/n
    end
    return rNew
end

#given two set of roots, we guess the new roots based on those
function smartGuess2ndorder(oldroots1, oldroots2,newN)
    n1,n2 = length(oldroots1), length(oldroots2)
    rRef1= linear_interpolation(range(1,n1),oldroots1)
    rRef2 = linear_interpolation(range(1,n2),oldroots2)
    rNew = zeros(newN)
    for i in range(1,newN)
        first= newN*rRef2(1+BigFloat((i-1)*(n2-1))/BigFloat(newN-1))/n2
        second = (newN-n2)* rRef2(1+BigFloat((i-1)*(n2-1))/BigFloat(newN-1) )
         - rRef1(1+BigFloat((i-1)*(n1-1))/BigFloat((newN-1)))/(n2-n1)
        rNew[i] = first+second

    end
   return rNew
end


#tab here is the list of Ns corresponding to number of roots
function better_guess(tab, root_list, newN; order=5)
    n=length(tab)
    new_x = vcat([1],[newN], [BigFloat(log(newN)^k) for k in 1:order])
    rfit= ones(BigFloat, (newN,))
    for i in range(1,newN)
        y=ones(BigFloat,(n,1))
        for j in range(1,n)
            l = length(root_list[tab[j]])
            p = linear_interpolation(range(1,l) , root_list[tab[j]][:,1], extrapolation_bc=Line() )
            y[j]=  p(BigFloat(i*tab[j])/BigFloat(newN))
        end

        X = hcat(ones(n),[x for x in tab], [BigFloat(log(x)^k) for x in tab,k in 1:order ])
        # print(size([ log(x)^k for x in tab,k in 1:order ]), size([x for x in tab]),size(ones(n)))
        
        rfit[i] = (transpose((X\y))*new_x)[1] 
    end
    return rfit
end

#We make prediction for the next root from the list of known roots labeled by tab
function better_guess2(tab, full_list, newN; order=5)
    n=length(tab)
    new_x = vcat([1],[newN], [BigFloat(log(newN)^k) for k in 1:order])
    rfit= ones(BigFloat, (newN,))
    for i in range(1,newN)
        y=ones(BigFloat,(n,1))
        for j in range(1,n)
            elem_length = length(full_list[tab[j]])
            l = Int(elem_length/2)
            p = linear_interpolation(range(1,l) , full_list[tab[j]][1:l,1], extrapolation_bc=Line() )
            y[j]=  p(BigFloat(i*tab[j])/BigFloat(newN))
        end

        X = hcat(ones(n),[x for x in tab], [BigFloat(log(x)^k) for x in tab,k in 1:order ])
        # print(size([ log(x)^k for x in tab,k in 1:order ]), size([x for x in tab]),size(ones(n)))
        
        rfit[i] = (transpose((X\y))*new_x)[1] 
    end
    return rfit
end


end

