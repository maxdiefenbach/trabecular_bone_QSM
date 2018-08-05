function deltaTE = calc_deltaTE(TE)

deltaTE = TE(2) - TE(1);
df = -220; cycTime = 1/abs(df);
if deltaTE > cycTime
    deltaTE = mod(deltaTE, cycTime);
else
    if deltaTE>3*cycTime/5,
        deltaTE = cycTime - deltaTE;
    end
end