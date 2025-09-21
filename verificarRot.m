function Rotacao = verificarRot(ordem, q)
    % Verifica o tipo de transformação baseado na ordem (x, y, z, Dx, Dy, Dz).
    % Retorna a matriz de rotação ou deslocamento correspondente.
    %
    % Entradas:
    %   ordem - Um caractere indicando o eixo ou tipo ('x', 'y', 'z', 'Dx', 'Dy', 'Dz')
    %   q     - O ângulo de rotação ou deslocamento
    %
    % Saída:
    %   Rotacao - Matriz de transformação homogênea (4x4)

    switch ordem
        case 'x'
            % Rotação em torno do eixo X
            Rotacao = rotacaoX(q);
        case 'y'
            % Rotação em torno do eixo Y
            Rotacao = rotacaoY(q);
        case 'z'
            % Rotação em torno do eixo Z
            Rotacao = rotacaoZ(q);
        case 'Dx'
            % Deslocamento variável ao longo de X
            Rotacao = Deslocamento([q, 0, 0]);
        case 'Dy'
            % Deslocamento variável ao longo de Y
            Rotacao = Deslocamento([0, q, 0]);
        case 'Dz'
            % Deslocamento variável ao longo de Z
            Rotacao = Deslocamento([0, 0, q]);
        otherwise
            % Exibir erro para ordens desconhecidas
            error('verificarRot: ordem desconhecida "%s". Use "x", "y", "z", "Dx", "Dy", ou "Dz".', ordem);
    end
end