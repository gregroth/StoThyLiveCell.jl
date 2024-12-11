
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

end