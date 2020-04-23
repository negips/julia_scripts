function specificGramian(A,B,tau,objTol)

    # Solves for the Gramian P:=int(exp(tA^T)Bexp(tA)dt)
    
    # Input:
    # A and B matrix, where A is anti-symmetric
    # tau: is the evaluation time
    # objTol: to set a limit for the number of terms to take in the Taylor expansion

    normB = norm(B);

    # Find the number of terms of the Taylor expansion to take. Consider only the leading order term
    ftol(x) = normB/factorial(Int32(round(x+1)))-objTol;
    err = 1000;
    N   = 0;

    # Can be done like this because the error is monotonically decreasing
    while err >= 0
        N +=1
        err = ftol(N);
    end

    println("To obtain an accuracy: ", objTol, " it requires: ",N," iterations");
    # Compute the necessary terms from the Taylor series
    # Do the first step and put the rest in a loop

    # for k=0
    P    = B*tau;
    Ck_1 = B;
    
    for k=1:N
        
        Ck = Ck_1*A-A*Ck_1;
        P  = P+tau^(k+1)*Ck/(factorial(k+1))
        Ck_1 = Ck;

    end

    return P

end



    
