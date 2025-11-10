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
                        -Dinamica.Izx, -Dinamica.Izy,  Dinamica.Izz ] .* Dinamica.zeros
    Dinamica.Rotacao    = Pcm.T(1:3,1:3);
    Dinamica.Posicao    = Pcm.Pcm;
    Dinamica.Velocidade = Pcm.Vcm;
    Dinamica.Omega      = Pcm.W;
    R   = Dinamica.Rotacao;
    w   = Dinamica.Omega(:);       % velocidade angular em base (3x1)
    v   = Dinamica.Velocidade(:);  % velocidade linear do CM em base (3x1)
    % Tensor de inércia do elo calculado no centro de massa e expresso na base
    Ic0 = R * Dinamica.Ip * R.';

    Dinamica.I0      = Ic0;
    Dinamica.EKrot   = (1/2) * (w.' * (Ic0 * w));
    Dinamica.EKtrans = (1/2) * Dinamica.m * (v.' * v);
    Dinamica.EK      = Dinamica.EKrot + Dinamica.EKtrans;

    Dinamica.EG         = -Dinamica.m * Dinamica.g * Dinamica.Posicao;
    Dinamica.L          = Dinamica.EK - Dinamica.EG;
end