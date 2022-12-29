function w = Wigner3j_Table(j1_max, j2_max, j3_max, m1, m2, m3)
    J1 = 169;
    J2 = 169;
    J3 = 59;
    % m1 = 0;
    % m2 = 0;
    % m3 = 0;

    fn = sprintf('Wigner3j_%d_%d_%d_%d_%d_%d.mat', J1, J2, J3, m1, m2, m3);
    w_table = load(fn).Expression1;

    tic
    w = zeros(J1+1, J2+1, J3+1);
    t = 0;
    for j1 = 0:J1
        for j2 = 0:J2
            for j3 = 0:J3
                if abs(j1-j2) <= j3 && j1+j2 >= j3
                    t = t + 1;
                    w(j1+1, j2+1, j3+1) = w_table(t);
                end
            end
        end
    end
    w = w(1:j1_max+1, 1:j2_max+1, 1:j3_max+1);
    
end