using Test

"""
Compare les résultats contenus dans deux tableaux différents sur les 14 instances
et retourne un vecteur des pires valeurs obtenues pour chacune de celles-ci.
"""
function find_worst(prim::Array,kruskal::Array)
    resultat = zeros(14)
    j = 1
    i = 1
    for i in 1:14
        max = -Inf
        for j in 1:4
            if prim[i,j] > max
                max = prim[i,j]
            end
            if kruskal[i,j] > max
                max = kruskal[i,j]
            end
        end
        resultat[i,j] = max
    end
    resultat
end

"""
Teste si toutes les valeurs obtenues sont à l'intérieur d'un multiple de 2
 de la valeur optimale pour chaque instance.
"""
function is_RSL_inside_UB(resultat::Array)
    rep = true
    if resultat[1] > 2 * 1610
        rep = false
    end
    if resultat[2] > 2 * 2020
        rep = false
    end
    if resultat[3] > 2 * 25395
        rep = false
    end
    # if resultat[4] > 2 * 1950
    #     rep = false
    # end
    if resultat[5] > 2 * 699
        rep = false
    end
    if resultat[6] > 2 * 937
        rep = false
    end
    if resultat[7] > 2 * 2085
        rep = false
    end
    if resultat[8] > 2 * 2707
        rep = false
    end
    if resultat[9] > 2 * 1272
        rep = false
    end
    if resultat[10] > 2 * 5046
        rep = false
    end
    if resultat[11] > 2 * 6942
        rep = false
    end
    if resultat[12] > 2 * 11461
        rep = false
    end
    if resultat[13] > 2 * 2763
        rep = false
    end
    if resultat[14] > 2 * 1273
        rep = false
    end
    rep
end

@test is_RSL_inside_UB(find_worst(resultats_prim_first_vs_last_vs_lightest_vs_heaviest(),resultats_kruskal_first_vs_last_vs_lightest_vs_heaviest()))
