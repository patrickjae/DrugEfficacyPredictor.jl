function copy(o::ModelConfiguration)
    c = typeof(o)()
    for f in fieldnames(typeof(o))
        setfield!(c, f, getfield(o, f))
    end
    c
end
