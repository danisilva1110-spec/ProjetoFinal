function Ry = rotacaoY(angulo_graus)
    % Esta função retorna a matriz de rotação em torno do eixo y em graus.
    c = cos(angulo_graus);
    s = sin(angulo_graus);
    Ry = [c, 0, s, 0; 
          0, 1, 0, 0; 
         -s, 0, c, 0; 
          0, 0, 0, 1];
end