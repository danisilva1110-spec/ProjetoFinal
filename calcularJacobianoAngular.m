function JacAng = calcularJacobianoAngular(robo, elo)
    % Verifica se ordem é uma célula de strings e converte para vetor de strings
    ordem = string(robo.ordem(1:elo));

    % Inicializa o Jacobiano Angular como uma matriz simbólica de zeros
    JacAng = sym(zeros(3, elo));

    % Percorre os elos para determinar os eixos de rotação
    for i = 1:elo
        % Rotação do frame anterior (i-1) expressa na base
        if i == 1
            Rprev = sym(eye(3));
        else
            Rprev = robo.T{i-1}(1:3,1:3);
        end
        switch ordem(i)
            case "x"
                JacAng(:,i) = Rprev(:,1);  % Primeiro eixo da rotação
            case "y"
                JacAng(:,i) = Rprev(:,2);  % Segundo eixo da rotação
            case "z"
                JacAng(:,i) = Rprev(:,3);  % Terceiro eixo da rotação
            case {"Dx", "Dy", "Dz"}  % Junta prismática -> não contribui para rotação
                JacAng(:,i) = [0; 0; 0];
            otherwise
                warning("Tipo de junta desconhecido: %s", ordem(i));
        end
    end
    
    % Simplificação simbólica do Jacobiano Angular
    JacAng = simplify(JacAng);
end
