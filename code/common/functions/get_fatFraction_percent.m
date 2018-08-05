function OutParams = get_fatFraction_percent(OutParams)
% OutParams = get_fatFraction_percent(OutParams)

    W = OutParams.water;
    F = OutParams.fat;
    
    whereWaterIsDominant = abs(W) > abs(F);
    
    FF = whereWaterIsDominant .* (1 - abs(W ./ (W + F))) + ...
         ~whereWaterIsDominant .* abs(F ./ (W + F));

    OutParams.fatFraction_percent = FF .* 100;

end