function mm_effectsize(base,candidate)
    (adjr2(candidate) - adjr2(base))/(1-adjr2(candidate))
end
