# query the dataframe with given parameters
function query(df,v,ℓ,R,K,δ,s,n,label)
    get_df = @chain df begin
        @rsubset :v == v
        @rsubset :ℓ == ℓ
        @rsubset :R == R
        @rsubset :K == K
        @rsubset :δ == δ
        @rsubset :s == s
        @rsubset :n == n
        @rsubset :label == label
    end
    return get_df
end

# simple unicode to latex converter
function unicode2latex(s)
    if s == "ℓ"
        return "interaction length \$R\$"
    elseif s == "R"
        return "system dimension \$R\$"
    elseif s == "K"
        return "intensity of counter-divergence effects \$K\$"
    elseif s == "δ"
        return "degree of noise \$\\eta\$"
    elseif s == "s"
        return "size of informed fish \$s\$"
    elseif s == "n"
        return "number of fish \$n\$"
    else
        return s
    end
end
