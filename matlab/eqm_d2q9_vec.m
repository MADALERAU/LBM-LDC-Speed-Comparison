function f_eq = eqm_d2q9_vec(rho, ksi, u, w)
    % rho: 1 x N_y x N_x
    % ksi: 9 x 2
    % u: 2 x N_y x N_x
    % w: 9 x 1
    c_s = 1/sqrt(3);

    f_eq = rho .* (1 + pagemtimes(ksi, u)/c_s^2 + pagemtimes(ksi, u).^2/(2*c_s^4) - sum(u.^2)/(2*c_s^2)).*w;

end