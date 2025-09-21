function dispCinematicaDireta(StructCinematicaDireta)

disp(['Transformação do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.T)

disp(['Jacobiano do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.JL)

disp(['Velocidade Linear do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.V)

disp(['Derivada do Jacobiano linear do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.dJL)

disp(['Aceleração Linear do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.A)

disp(['Jacobiano Angular do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.JA)

disp(['Velocidade Angular do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.W)

disp(['Derivada do Jacobiano Angular do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.dJA)

disp(['Aceleração angular do ponto ' StructCinematicaDireta.nome ':']);
pretty(StructCinematicaDireta.alfa)

end