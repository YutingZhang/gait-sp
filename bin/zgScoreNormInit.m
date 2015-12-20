function [ NormCallBack RSaN ] = zgScoreNormInit( RSa, PreNormFunName )

NormCallBack = @zgAssign;
switch PreNormFunName
    case 'none'
        RSaN = RSa;
    case 'exp'
        RSaN = exp(RSa);
    case 'neg'
        RSaN = -RSa;
        NormCallBack = @uminus;
    otherwise
        error('Unrecognized param.PreNormFunName');
end

end


function A = zgAssign(A)
end
