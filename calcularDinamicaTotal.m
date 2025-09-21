function DinamicaTotal = calcularDinamicaTotal(Pcm,zeros)
    DinamicaTotal         = struct();
    DinamicaTotal.n       = length(Pcm);
    DinamicaTotal.zeros   = zeros;
    DinamicaTotal.dinCell = cell(1, DinamicaTotal.n);
    DinamicaTotal.EKTotal = 0;
    DinamicaTotal.EGTotal = 0;
    DinamicaTotal.LTotal  = 0;
    for i=1:DinamicaTotal.n
        DinamicaTotal.dinCell{i} = calcularDinamicaDosPontos(Pcm{i},DinamicaTotal.zeros{i});
        DinamicaTotal.EKTotal    = DinamicaTotal.EKTotal + DinamicaTotal.dinCell{i}.EK;
        DinamicaTotal.EGTotal    = DinamicaTotal.EGTotal + DinamicaTotal.dinCell{i}.EG;
        DinamicaTotal.LTotal     = DinamicaTotal.LTotal  + DinamicaTotal.dinCell{i}.L;
    end
    
    DinamicaTotal.Torque = calcularTorque(Pcm, DinamicaTotal);

end