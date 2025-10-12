function Dinamica = calcularDinamica(Pcm)
    syms g

    n = Pcm{end}.n;
    syms Massa [1,n]
    syms Ixx   [n,1]
    syms Ixy   [n,1]
    syms Ixz   [n,1]
    syms Iyx   [n,1]
    syms Iyy   [n,1]
    syms Iyz   [n,1]
    syms Izx   [n,1]
    syms Izy   [n,1]
    syms Izz   [n,1]

    qs  = Pcm{end}.juntas(1,:);
    dqs = Pcm{end}.juntas(2,:);

    I = cell(1,n);
    for i = 1:n
        I{i} = [ Ixx(i), -Ixy(i), -Ixz(i);
                -Iyx(i),  Iyy(i), -Iyz(i);
                -Izx(i), -Izy(i),  Izz(i)];
    end

    Dinamica           = struct();
    Dinamica.G         = sym(zeros(n,1));
    Dinamica.M         = sym(zeros(n,n));
    Dinamica.H1        = sym(zeros(n,1));
    Dinamica.H2        = sym(zeros(n,1));
    Dinamica.Ec        = sym(zeros(1,n));
    Dinamica.Ep        = sym(zeros(1,n));
    Dinamica.EcT       = sym(0);
    Dinamica.EpT       = sym(0);
    Dinamica.gravidade = [0 0 g];

    for i = 1:n
        Raux           = Pcm{i}.T(1:3,1:3);
        Iaux           = Raux * I{i} * transpose(Raux);
        Dinamica.Ep(i) = Massa(i) * Dinamica.gravidade * (Pcm{i}.P);
        Dinamica.Ec(i) = (Massa(i)/2) * transpose(Pcm{i}.V) * Pcm{i}.V + ...
                         (1/2) * transpose(Pcm{i}.W) * Iaux * (Pcm{i}.W);
    end

    Dinamica.EcT = sum(Dinamica.Ec);
    Dinamica.EpT = sum(Dinamica.Ep);

    Dinamica.Lagrangeano = Dinamica.EcT - Dinamica.EpT;

    for i = 1:n
        Dinamica.G(i,1)  =  diff(Dinamica.EpT,qs(i));
        Dinamica.H2(i,1) = -diff(Dinamica.EcT,qs(i));
        for j = 1:n
            Dinamica.M(i,j) = diff(diff(Dinamica.EcT,dqs(i)),dqs(j));
            Dinamica.H1(i,1) = Dinamica.H1(i,1) + diff(diff(Dinamica.EcT,dqs(i)),qs(j)) * dqs(j);
        end
    end

    Dinamica.H = Dinamica.H1 + Dinamica.H2;
end
