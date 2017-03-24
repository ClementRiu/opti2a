
function [F] = critere(q_c)

    qbq = q0 + B * q_c;
    F = (qbq' * (r .* qbq .* abs(qbq)) / 3 + (pr' * (Ar*qbq)));

endfunction

function [G] = gradient(q_c)

    qbq = q0 + B * q_c;

    G = B' * ( (r .* abs(q0 + B * q_c) .* (q0 + B * q_c)) + Ar' * pr);

endfunction

function [H] = hessien(q_c)

    qbq = q0 + B * q_c;

    H = 2 * B' * (eye(q0 * q0') .* (r * abs(qbq)')) * B;

endfunction

function [F, G, ind]=OraclePG(q_c,ind)
    F = 0;
    G = 0;

    if ind == 2  | ind==4 then
        F = critere(q_c);
    end

    if ind == 3 | ind==4 then
        G = gradient(q_c);
    end

endfunction

function [F, G, H, ind]=OraclePH(q_c,ind)
    F = 0;
    G = 0;
    H = 0;


    if ind == 2  | ind==4 | ind == 7 then
        F = critere(q_c);
    end

    if ind == 3 | ind==4 | ind ==6 | ind == 7 then
        G = gradient(q_c);
    end

    if ind == 5 | ind == 6 | ind == 7 then
        H = hessien(q_c)
    end

endfunction
