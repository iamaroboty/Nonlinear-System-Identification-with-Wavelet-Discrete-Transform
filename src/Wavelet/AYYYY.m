    % Decompose approximation and detail coefficient
    low_d = 1:4;
    d = 125;
    lf = length(low_d);
    LL = [d; zeros(level,1)];
    for i= 1:level
        LL = [floor((LL(1)+lf-1)/2); LL(1:end-1)];
    end
    LL = [LL(1); LL]';

    % Z = [cA||cD||...||cD]
    zA = Z(1:L(1));
    for i=1:level
        zD{i} = Z(L(i)+1:sum(L(i),L(i+1)));
    end