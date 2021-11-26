function AtomBase.diagonal(op::Operator{T,S}) where {T<:Quantity,S<:AtomState}
    units = unit(T)
    opu = Operator(ustrip.(op.c), op.ks, op.bs)
    vals, vecs = AtomBase.diagonal(opu)
    return (; values = vals * units, vectors = vecs)
end