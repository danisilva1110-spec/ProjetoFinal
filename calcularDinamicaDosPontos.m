function Dinamica = calcularDinamicaDosPontos(Pcm,zeros)
    Dinamica         = struct();
    Dinamica.zeros   = zeros;
    Dinamica.elo     = num2str(size(Pcm.juntas, 2));
    Dinamica.gmod    = sym(['g'   Dinamica.elo]);    
    Dinamica.m       = sym(['m'   Dinamica.elo]);
    Dinamica.Ixx     = sym(['Ixx' Dinamica.elo]);
    Dinamica.Ixy     = sym(['Ixy' Dinamica.elo]);
    Dinamica.Ixz     = sym(['Ixz' Dinamica.elo]);
    Dinamica.Iyx     = sym(['Iyx' Dinamica.elo]);
    Dinamica.Iyy     = sym(['Iyy' Dinamica.elo]);
    Dinamica.Iyz     = sym(['Iyz' Dinamica.elo]);
    Dinamica.Izx     = sym(['Izx' Dinamica.elo]);
    Dinamica.Izy     = sym(['Izy' Dinamica.elo]);
    Dinamica.Izz     = sym(['Izz' Dinamica.elo]);
    Dinamica.g       = [            0,             0, -Dinamica.gmod];
    Dinamica.Ip      = [ Dinamica.Ixx, -Dinamica.Ixy, -Dinamica.Ixz;
                        -Dinamica.Iyx,  Dinamica.Iyy, -Dinamica.Iyz;
                        -Dinamica.Izx, -Dinamica.Izy,  Dinamica.Izz ] .* Dinamica.zeros;
    Dinamica.Rotacao    = Pcm.T(1:3,1:3);
    Dinamica.Posicao    = Pcm.Pcm;
    Dinamica.Velocidade = Pcm.Vcm;
    Dinamica.Omega      = Pcm.W;
    Dinamica.I0         = Dinamica.Rotacao*Dinamica.Ip*transpose(Dinamica.Rotacao);
    Dinamica.EK         = (1/2)*Dinamica.m*transpose(Dinamica.Velocidade)*Dinamica.Velocidade + (1/2)*transpose(Dinamica.Omega)*(Dinamica.I0*Dinamica.Omega);
    Dinamica.EG         = -Dinamica.m * Dinamica.g * Dinamica.Posicao;
    Dinamica.L          = Dinamica.EK - Dinamica.EG;
end