function linPhaseMat = get_linPhaseMat(alpha, ht, wd)
linPhaseMat = exp(i.*repmat(((-ht/2+1):ht/2)'*alpha, [1 wd]));
return;