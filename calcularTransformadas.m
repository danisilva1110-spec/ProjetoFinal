function Transformadas = calcularTransformadas(ordem, excentricidades)

    %% Criação das variáveis simbólicas
    Transformadas   = struct();
    Transformadas.n = length(ordem);
    syms Q   [1, Transformadas.n];  % Posições articulares
    syms dQ  [1, Transformadas.n];  % Velocidades articulares
    syms d2Q [1, Transformadas.n];  % Acelerações articulares
    syms Lx  [1, Transformadas.n];  % Comprimentos ao longo de X
    syms Ly  [1, Transformadas.n];  % Comprimentos ao longo de Y
    syms Lz  [1, Transformadas.n];  % Comprimentos ao longo de Z
    
    Transformadas.L   = [transpose(Lx), transpose(Ly), transpose(Lz)];
    Transformadas.L   = Transformadas.L .* excentricidades; 
    Transformadas.ordem = ordem;

    %% Cálculo das transformadas   
    Taux   = eye(4,4);  % Transformação inicial
    Transformadas.T   = cell(1, Transformadas.n);  % Lista de transformações

    for i = 1:Transformadas.n
        % Verificar se é uma rotação ou deslocamento
        if contains(ordem{i}, 'D')
            % Comprimento variável (parametrizado por Q(i))
            deslocamento = [0, 0, 0];
            if ordem{i} == "Dx"
                deslocamento = [Q(i), 0, 0];  % Deslocamento em X
            elseif ordem{i} == "Dy"
                deslocamento = [0, Q(i), 0];  % Deslocamento em Y
            elseif ordem{i} == "Dz"
                deslocamento = [0, 0, Q(i)];  % Deslocamento em Z
            end
            Taux = Taux * Deslocamento(deslocamento);
        else
            % Caso seja uma rotação com comprimento fixo
            Taux = Taux * verificarRot(ordem{i}, Q(i)) * Deslocamento(Transformadas.L(i, :));
        end

        % Salvar a transformação acumulada
        Transformadas.T{i} = Taux;
    end

    % Juntas: posição, velocidade e aceleração
    Transformadas.juntas = [Q; dQ; d2Q];
end