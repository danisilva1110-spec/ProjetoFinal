function cinTotal = calcularCinematicaTotal(robo)
    cinTotal = cell(1, robo.n);
    for i=1:robo.n
        cinTotal{i}=cinematicaDireta(robo, i);
    end
end