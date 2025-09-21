 function dJ = DerivadaJacobiano(Jacobiano,Q,derivadaQ)
    [linha, coluna] = size(Jacobiano);
    dJ = sym('zero', [linha, coluna]);;
    for n=1:coluna
        dJ(:,n)=jacobian(transpose(Jacobiano(:,n)),Q)*transpose((derivadaQ));
    end
end