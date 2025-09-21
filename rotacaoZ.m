function Rz = rotacaoZ(angulo_graus)
    % Esta função retorna a matriz de rotação em torno do eixo z em graus.
    c = cos(angulo_graus);
    s = sin(angulo_graus);
    Rz = [c,-s, 0, 0; 
          s, c, 0, 0; 
          0, 0, 1, 0; 
          0, 0, 0, 1];
end