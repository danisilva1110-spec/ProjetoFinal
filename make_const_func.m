function make_const_func(funcName, value)
% Cria um .m do tipo:
%   function X = funcName()
%       X = [ ...numeros... ];
%   end
%
% Ex.: make_const_func('funcM', Dinamica.Torque.M);

    if ~isnumeric(value) && ~islogical(value)
        error('Apenas arrays numéricos/lógicos.');
    end
    fid = fopen([funcName '.m'],'w');
    assert(fid>0, 'Não consegui criar o arquivo.');

    fprintf(fid,'function X = %s()\n', funcName);
    % 16 dígitos para preservar precisão
    fprintf(fid,'X = %s;\n', mat2str(value, 16));
    fprintf(fid,'end\n');
    fclose(fid);
    fprintf('Gerado: %s.m\n', funcName);
end
