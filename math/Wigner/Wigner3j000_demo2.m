clear all

J1 = 2;
J2 = 2;
J3 = 2;

tic
t = 0;
for j1 = 0:J1
    for j2 = 0:J2
        for j3 = 0:J3
            if abs(j1-j2) <= j3 && j1+j2 >= j3
                t = t + 1;
            end
        end
    end
end
toc

(J1+1)*(J2+1)*(J3+1)
t
