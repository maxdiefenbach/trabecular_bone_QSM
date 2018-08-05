function dTE_s = get_dTE_s(ImDataParams)

    TE_s = ImDataParams.TE_s;
    dTE_s = TE_s(2) - TE_s(1);

end