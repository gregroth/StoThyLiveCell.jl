
module utiles

export mergeparameter, mergeparameter_base

#merge the free and fixed parameters

function mergeparameter(fixedparameter, freeparameter, freeparameteridx)
    if fixedparameter[1] == -1.
        return freeparameter
    else
        nbparam = length(fixedparameter) + length(freeparameter)
        parameter = Vector{Float64}(undef,nbparam)
        kfree = 1
        kfixed = 1
        @inbounds for i =1 : nbparam
            if sum((i .== freeparameteridx)) == 1
                parameter[i] = freeparameter[kfree]
                kfree += 1
            else
                parameter[i] = fixedparameter[kfixed]
                kfixed += 1
            end
        end
        return parameter
    end
end    

function mergeparameter_base(fixedparameter, freeparameter::AbstractVector{T}, freeparameteridx) where T
    if fixedparameter[1] == -1
        return freeparameter
    else
        nbparam = length(fixedparameter) + length(freeparameter)
        parameter = Vector{T}(undef,nbparam)
        kfree = 1
        kfixed = 1
        for i =1 : nbparam
            if sum((i .== freeparameteridx)) == 1
                parameter[i] = freeparameter[kfree]
                kfree += 1
            else
                parameter[i] = fixedparameter[kfixed]
                kfixed += 1
            end
        end
        return parameter
    end
end    

function convolution(dist::Vector{Float64})
    conv = zeros(Float64, 2*length(dist) - 1)
    for i in eachindex(dist)
        for j in eachindex(dist)
            conv[i + j - 1] += dist[i] * dist[j]
        end
    end
    return conv
end

end