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
    r   = Dinamica.Posicao(:);     % CM em base (3x1)
    w   = Dinamica.Omega(:);       % velocidade angular em base (3x1)
    Ic0 = R * Dinamica.Ip * R.';   % tensor no CM expresso na base
    Ipar = Dinamica.m * ( (r.'*r)*eye(3) - (r*r.') );  % eixo paralelo (matricial)
    I0_about_base = Ic0 + Ipar;
    
    Dinamica.I0 = I0_about_base;
    Dinamica.EK = (1/2) * (w.' * (I0_about_base * w));   % <-- sem termo translacional

    Dinamica.EG         = -Dinamica.m * Dinamica.g * Dinamica.Posicao;
    Dinamica.L          = Dinamica.EK - Dinamica.EG;
end