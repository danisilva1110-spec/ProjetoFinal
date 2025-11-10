function Transformadas = calcularTransformadas(ordem, excentricidades)
    %% Criação das variáveis simbólicas
    Transformadas           = struct();
    Transformadas.n         = length(ordem);
    Transformadas.excentricidades = excentricidades;

    syms Q   [1, Transformadas.n];  % Posições articulares
    syms dQ  [1, Transformadas.n];  % Velocidades articulares
    syms d2Q [1, Transformadas.n];  % Acelerações articulares
    syms Lx  [1, Transformadas.n];  % Comprimentos ao longo de X (fixos do elo)
    syms Ly  [1, Transformadas.n];  % Comprimentos ao longo de Y
    syms Lz  [1, Transformadas.n];  % Comprimentos ao longo de Z

    % Centro de massa local do elo i (em relação ao frame do elo i)
    % Observação: Cx(i), Cy(i), Cz(i) podem ser funções simbólicas de Q(i), se desejar.
    syms Cx  [1, Transformadas.n];
    syms Cy  [1, Transformadas.n];
    syms Cz  [1, Transformadas.n];

    Transformadas.L     = [transpose(Lx), transpose(Ly), transpose(Lz)];
    Transformadas.L     = Transformadas.L .* Transformadas.excentricidades;
    Transformadas.ordem = ordem;

    Transformadas.CM    = [transpose(Cx), transpose(Cy), transpose(Cz)]; % n x 3

    %% Cálculo das transformadas
    Taux                = eye(4,4);                 % Transformação acumulada (base -> frame atual)
    Transformadas.T     = cell(1, Transformadas.n); % Base -> ponta do elo i
    Transformadas.Tcm   = cell(1, Transformadas.n); % Base -> CM do elo i
    Transformadas.r_cm  = cell(1, Transformadas.n); % Vetor posição do CM do elo i (em base)

    for i = 1:Transformadas.n
        oi = ordem{i};

        if contains(oi, 'D')  % Junta prismática (deslocamento variável)
            deslocamento = [0, 0, 0];
            if oi == "Dx"
                deslocamento = [Q(i), 0, 0];
            elseif oi == "Dy"
                deslocamento = [0, Q(i), 0];
            elseif oi == "Dz"
                deslocamento = [0, 0, Q(i)];
            end

            % --- Transformação até o CM do elo i ---
            % CM após o deslocamento variável + offset local de CM
            Tcm_i = Taux * Deslocamento(deslocamento + Transformadas.CM(i,:));
            Transformadas.Tcm{i}  = Tcm_i;
            Transformadas.r_cm{i} = Tcm_i(1:3,4);

            % --- Atualiza a transformação até a "ponta" do elo i ---
            % Se existir um comprimento fixo L(i,:) do elo após a junta prismática
            Taux = Taux * Deslocamento(deslocamento) * Deslocamento(Transformadas.L(i,:));

        else % Junta rotacional (Rx, Ry, Rz, etc.) com comprimento fixo L(i,:)
            Troti = verificarRot(oi, Q(i));

            % --- Transformação até o CM do elo i ---
            % CM após a rotação da junta i, antes do deslocamento total L(i,:)
            Tcm_i = Taux * Troti * Deslocamento(Transformadas.CM(i,:));
            Transformadas.Tcm{i}  = Tcm_i;
            Transformadas.r_cm{i} = Tcm_i(1:3,4);

            % --- Atualiza a transformação até a "ponta" do elo i ---
            Taux = Taux * Troti * Deslocamento(Transformadas.L(i,:));
        end

        % Salva a transformação acumulada até a ponta do elo i
        Transformadas.T{i} = Taux;
    end

    % Juntas: posição, velocidade e aceleração
    Transformadas.juntas = [Q; dQ; d2Q];
end
