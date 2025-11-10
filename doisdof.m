matrizExcentricidades = [1 0 0;   % 1) p1: pedestal -> eixo 2  (ex.: d1 ao longo de Z)
                         1 0 0];   % 2) p2: elo 1 -> elo 2     (ex.: a2 ao longo de X)

ordem = {'z','z'};

matrizNaBaseInercial = calcularTransformadas(ordem, matrizExcentricidades);

pontosDeJuntas = calcularCinematicaTotal(matrizNaBaseInercial);

dispCinematicaDireta(pontosDeJuntas{2});

Izero1 = [0 0 0;
          0 0 0;
          0 0 1];
Izero2 = [0 0 0;
          0 0 0;
          0 0 1];

Izeros = {Izero1, Izero2};

Dinamica = calcularDinamicaTotal(pontosDeJuntas, Izeros);