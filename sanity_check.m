% ========== 1) Monta robô antropomórfico ==========
ordem = {'z','y','y','z','y','z'};
matrizExcentricidades = [ 0 0 1;   % elo 1: Z
                          1 0 0;   % elo 2: X
                          1 0 0;   % elo 3: X
                          0 0 1;   % elo 4: Z
                          0 0 1;   % elo 5: Z
                          0 0 1 ]; % elo 6: Z

robo = calcularTransformadas(ordem, matrizExcentricidades);

% Símbolos (forma compatível)
Q   = sym('Q',   [1 6], 'real');
dQ  = sym('dQ',  [1 6], 'real');
d2Q = sym('d2Q', [1 6], 'real');
Lx  = sym('Lx',  [1 6], 'real');
Ly  = sym('Ly',  [1 6], 'real');
Lz  = sym('Lz',  [1 6], 'real');

% Centros de massa simbólicos (um por elo)
Lxcm = sym('Lxcm', [1 6], 'real');
Lycm = sym('Lycm', [1 6], 'real');
Lzcm = sym('Lzcm', [1 6], 'real');

% Comprimentos numéricos de exemplo
Lvals = [ 0    0    0.10;  % elo 1
          0.20 0    0.00;  % elo 2
          0.15 0    0.00;  % elo 3
          0    0    0.10;  % elo 4
          0    0    0.08;  % elo 5
          0    0    0.05]; % elo 6

% Pose/vel/acc de teste
Qv   = [  0.2  -0.3   0.4  -0.2   0.1  -0.1 ];
dQv  = [  0.5   0.3   0.0   0.2  -0.1   0.05];
d2Qv = [  0.2  -0.1   0.1   0.0   0.0   0.0 ];

% Centros de massa numéricos (exemplo simples)
LcmVals = repmat([0.05 0 0], 6, 1);  % ajuste conforme seus dados

% Listas de substituição (inclui Lxcm/Lycm/Lzcm!)
subsList = [ num2cell(Q),  num2cell(dQ),  num2cell(d2Q), ...
             num2cell(Lx), num2cell(Ly), num2cell(Lz), ...
             num2cell(Lxcm), num2cell(Lycm), num2cell(Lzcm) ];

valsList = [ num2cell(Qv),  num2cell(dQv),  num2cell(d2Qv), ...
             num2cell(Lvals(:,1).'), num2cell(Lvals(:,2).'), num2cell(Lvals(:,3).'), ...
             num2cell(LcmVals(:,1).'), num2cell(LcmVals(:,2).'), num2cell(LcmVals(:,3).') ];

% ========== 2) Cinemática direta até cada elo ==========
cin = calcularCinematicaTotal(robo);   % cell 1x6 (usa cinematicaDireta)

% ========== 3) Checagens JLcm e Pcm ==========
for i = 1:6
    c = cin{i};

    % (A) JLcm vs jacobian(Pcm,q)
    JL_from_def = jacobian((c.Pcm).', c.q);
    JL_num  = double(subs(c.JLcm, subsList, valsList));
    JL2_num = double(subs(JL_from_def, subsList, valsList));
    errJL = norm(JL_num - JL2_num, 'fro');

    % (B) Pcm = T{elo} * Deslocamento(Lcm)
    Tcur_num  = double(subs(robo.T{i}, subsList, valsList));
    Lcm_i     = LcmVals(i,:).';
    Pcm_fromT = Tcur_num * [eye(3) Lcm_i; 0 0 0 1];
    Pcm_fromT = Pcm_fromT(1:3,4);
    Pcm_here  = double(subs(c.Pcm, subsList, valsList));
    errPcm = norm(Pcm_here - Pcm_fromT);

    fprintf('Elo %d: ||JL-Jacobian(Pcm)|| = %.3e   ||Pcm - T*Lcm|| = %.3e\n', i, errJL, errPcm);
end

% ========== 4) Dinâmica: M simétrica e G coerente ==========
mval = ones(1,6);
Ival = 1e-2*ones(6,1);  % inércia principal fictícia
Izeros = repmat({ones(3)}, 1, 6);

Dcm = calcularDinamicaTotal(cin, Izeros);

subsD = subsList; valsD = valsList;
for i = 1:6
    subsD = [subsD, {sym(sprintf('m%d',i))}];      valsD = [valsD, {mval(i)}];
    subsD = [subsD, {sym(sprintf('g%d',i))}];      valsD = [valsD, {9.81}];
    subsD = [subsD, {sym(sprintf('Ixx%d',i)), sym(sprintf('Iyy%d',i)), sym(sprintf('Izz%d',i)), ...
                     sym(sprintf('Ixy%d',i)), sym(sprintf('Ixz%d',i)), sym(sprintf('Iyz%d',i)), ...
                     sym(sprintf('Iyx%d',i)), sym(sprintf('Izx%d',i)), sym(sprintf('Izy%d',i))}];
    valsD = [valsD, {Ival(i), Ival(i), Ival(i), 0, 0, 0, 0, 0, 0}];
end

Mnum = double(subs(Dcm.Torque.M, subsD, valsD));
Gnum = double(subs(Dcm.Torque.G, subsD, valsD));

fprintf('\n||M - M^T||_F = %.3e\n', norm(Mnum - Mnum.', 'fro'));
fprintf('Norma de G: %.3e\n', norm(Gnum));
detM = det(Mnum);
fprintf('det(M) = %.3e\n', detM);

if norm(Mnum - Mnum.', 'fro') < 1e-8
    disp('✅ M simétrica (ok)');
else
    warning('⚠️  M não está exatamente simétrica — verifique as expressões.');
end

if abs(detM) > 1e-8 && isfinite(detM)
    disp('✅ determinante de M diferente de zero (ok)');
else
    warning('⚠️  determinante de M está próximo de zero ou não numérico.');
end
