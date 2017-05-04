function [F] = critere(lagrange)
       
    z = (Ar' * pr + Ad' * lagrange) ./ r;
    
    sgn = ((z <= 0) - 0.5) * 2

    qlamb = sqrt(abs(z)) .* sgn;
    F = (qlamb' * (r .* qlamb .* abs(qlamb)) / 3 + (pr' * (Ar * qlamb))) + (Ad * qlamb - fd)' * lagrange;
    F = -F;

endfunction

function [G] = gradient(lagrange)

    z = (Ar' * pr + Ad' * lagrange) ./ r;
    
    sgn = ((z <= 0) - 0.5) * 2

    qlamb = sqrt(abs(z)) .* sgn;

    G = Ad * qlamb - fd;
    G = -G;

endfunction

function [H] = hessien(lagrange)
    
    
    y = (Ar' * pr + Ad' * lagrange);
    
    z = y ./ r;
    
    sgny = ((y >= 0) - 0.5) * 2;
    sgn = ((z <= 0) - 0.5) * 2;
    sgnr = ((r >= 0) - 0.5) * 2; 
 
    H = Ad * diag((sign(-z) ./ sqrt(abs(r))) .* (sign(y) ./ sqrt(abs(y)))) * Ad' / 2;

    H = -H;

endfunction

function [compteurF, G, ind]=OraclePGL(lagrange, ind, compteur)
    compteur = compteur + 1;
    
    F = 0;
    G = 0;

    if ind == 2  | ind==4 then
        F = critere(lagrange);
    end

    if ind == 3 | ind==4 then
        G = gradient(lagrange);
    end

endfunction

function [compteur, F, G, H, ind]=OraclePHL(lagrange, ind, compteur)
    compteur = compteur + 1;
    
    F = 0;
    G = 0;
    H = 0;


    if ind == 2  | ind==4 | ind == 7 then
        F = critere(lagrange);
    end

    if ind == 3 | ind==4 | ind ==6 | ind == 7 then
        G = gradient(lagrange);
    end

    if ind == 5 | ind == 6 | ind == 7 then
        H = hessien(lagrange)
    end

endfunction
