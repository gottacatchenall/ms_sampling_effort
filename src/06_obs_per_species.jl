using EcologicalNetworks
using Distributions


S = 50


for relabd in 0:0.01:0.1
    for numobs in [10, 50, 100, 250, 500, 1000]
        for rep in 1:nreps
            
            # use a gradient of relative abundance as simlate
            # as probability of observing as bernoulli with p = relabd.

            observedspecies = zeros(S)
            while thisspeciesobs < numobs
                observedspecies = rand(Categorical(relabd))   
                if observedspecies == thisspecies
                    thisspeciesobs += 1
                end
            end  
            fn, tp = 0, 0
            for i in 1:S, j in 1:S
                tp += A[i,j] && i in observations && j in observations
                fn += A[i,j] && i in observations || j in observations
            end

            thisfnr = fn/(fn+tp)
        end
    end
end
