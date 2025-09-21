function D = Deslocamento(VetorDeslocamento)
    D = [ 1, 0, 0, VetorDeslocamento(1);
          0, 1, 0, VetorDeslocamento(2);
          0, 0, 1, VetorDeslocamento(3);
          0, 0, 0, 1                    ];
end