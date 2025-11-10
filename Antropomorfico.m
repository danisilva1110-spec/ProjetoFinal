% === Excentricidades (quais eixos têm comprimento fixo por junta) ===
% Lx, Ly, Lz serão criados simbolicamente em calcularTransformadas.
% Aqui só marcamos em qual eixo cada junta tem deslocamento fixo (0/1).

matrizExcentricidades = [0 0 1;   % 1) p1: pedestal -> eixo 2  (ex.: d1 ao longo de Z)
                         1 0 0;   % 2) p2: elo 1 -> elo 2     (ex.: a2 ao longo de X)
                         1 0 0;   % 3) p3: elo 2 -> elo 3     (ex.: a3 ao longo de X)
                         1 0 0;   % 4) p4: elo 3 -> wrist center (absorve t_{3->wc} no T4)
                         0 0 0;   % 5) p5: punho esférico (rotação pura no mesmo ponto)
                         0 0 1];  % 6) p6: wrist center -> TCP (absorve ferramenta em Z)

% === Ordem das juntas (apenas 6 transformações) ===
% Cada T_i = Rotacao_eixo(q_i) * Trans(p_i)
ordem = {'z','y','y','z','y','z'};

% Calcular as transformadas homogêneas para a nova configuração
matrizNaBaseInercial = calcularTransformadas(ordem, matrizExcentricidades);

% Calcular os pontos das juntas a partir das transformações
pontosDeJuntas = calcularCinematicaTotal(matrizNaBaseInercial);

% Exibir o ponto final (cinemática direta)
dispCinematicaDireta(pontosDeJuntas{6});

% Matrizes de inércia para as três rotações em Z
Izero1 = [1 0 0;
          0 1 0;
          0 0 1];
Izero2 = [1 0 0;
          0 1 0;
          0 0 1];
Izero3 = [1 0 0;
          0 1 0;
          0 0 1];
Izero4 = [1 0 0;
          0 1 0;
          0 0 1];
Izero5 = [1 0 0;
          0 1 0;
          0 0 1];
Izero6 = [1 0 0;
          0 1 0;
          0 0 1];

% Lista completa de inércias
Izeros = {Izero1, Izero2, Izero3, Izero4, Izero5, Izero6};

% Calcular a dinâmica total do robô
Dinamica = calcularDinamicaTotal(pontosDeJuntas, Izeros);