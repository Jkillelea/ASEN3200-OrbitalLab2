function v = gibbs(r1, r2, r3, mu, use_vector_num)
    n_r1 = norm(r1);
    n_r2 = norm(r2);
    n_r3 = norm(r3);

    N = n_r1*cross(r2, r3) + n_r2*cross(r3, r1) + n_r3*cross(r1, r2);
    D = cross(r1, r2) + cross(r2, r3) + cross(r3, r1);
    S = r1*(n_r2 - n_r3) + r2*(n_r3 - n_r1) + r3*(n_r1 - n_r2);

    switch use_vector_num
    case 1
        r = r1;
        n_r = n_r1;
    case 2
        r = r2;
        n_r = n_r2;
    case 3
        r = r3;
        n_r = n_r3;
    end

    v = sqrt(mu/(norm(N)*norm(D)))*(cross(D, r)/n_r + S);
end

