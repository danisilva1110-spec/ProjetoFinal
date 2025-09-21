function Torque = calcularTorque(Pcm, Dcm)
    Torque        = struct();
    Torque.n      = length(Pcm);
    Torque.M      = sym(zeros(Torque.n, Torque.n));
    Torque.qs     = Pcm{Torque.n}.juntas(1,:);
    Torque.dqs    = Pcm{Torque.n}.juntas(2,:);
    Torque.d2qs   = Pcm{Torque.n}.juntas(3,:);
    for i = 1:Torque.n
        for j = 1:Torque.n
            Torque.M(i,j) = simplify(diff(diff(Dcm.EKTotal, Torque.dqs(i)),Torque.dqs(j)), 'Steps', 5);
        end
    end
    
    Torque.C = sym(zeros(Torque.n, 1));
    Torque.H = sym(zeros(Torque.n, 1));

    for i = 1:Torque.n
        C_sum = 0;
        H_sum = 0;
        for j = 1:Torque.n
            for k = 1:Torque.n
                Ci_jk = (1/2) * ( diff(Torque.M(i,j), Torque.qs(k)) + diff(Torque.M(i,k), Torque.qs(j)) - diff(Torque.M(j,k), Torque.qs(i)) );
                if j == k
                    H_sum = H_sum + Ci_jk * Torque.dqs(j)^2;
                else
                    C_sum = C_sum + Ci_jk * Torque.dqs(j) * Torque.dqs(k);
                end                
            end
        end
        Torque.C(i) = simplify(C_sum, 'Steps', 5);
        Torque.H(i) = simplify(H_sum, 'Steps', 5);        
    end

    Torque.G = sym(zeros(Torque.n,1));

    for i = 1:Torque.n
        Torque.G(i) = simplify(diff(Dcm.EGTotal, Torque.qs(i)), 'Steps', 5);
    end
    
    Torque.Tau = Torque.M * transpose(Torque.d2qs) + Torque.C + Torque.H + Torque.G;
   
end