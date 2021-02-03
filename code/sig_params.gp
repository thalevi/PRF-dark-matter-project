

/* How to use this script: start PARI-GP with the 'gp' command, type '\r sig_params.gp', then call find_params_mod2mod3_L1() */
/* If you don't have gp installed, on debian/ubuntu the package name is 'pari-gp' */

/* Base 2 logarithm */
log2(x)=
{
    return(log(x)/log(2));
}

epsilon(M,n,tao)=
{
    local(i, z, maxz);
    
    maxz = 10^(-100);
    for(i = M-tao, M, 
       z = binomial(i, M-tao)/(binomial(M, M-tao)*n^(i-M+tao));
        if(z > maxz, 
            maxz = z;
         );
    );
    
    \\print("Soundness = ", maxz);
    \\print("in bits: ", log(1/maxz)/log(2));
    return(log(1/maxz)/log(2));
}


get_size_mod2mod3(M, N, t)= 
{
    /* Note: t stands for tau, in the write-up's notation */

    local(bitlen, KB, n, k, m, seed, hash, salt, aux, msgs, inputs, treeN, treeM );
  
    /* OWF parameters */
    n = 128;        /* Input size, bits */
    m = 453;        /* Intermediate size, denoted m in the paper, bits */
    k = 81;         /* Output size, number of Z_3 elements */
    ell3 = 1.58;    /* Bits req'd to represent a Z_3 element */

    secpar = 128;       /* Security level (classical; quantum is secpar/2) */
    seed = secpar;      /* Size of a random seed */
    hash = 2*secpar;    /* Digest length for commitments */     
    salt = 256;         /* Bitlength of salt, always 256 */

    /* |MPC| for the Original (2,3)-OWF in Itai's writeup */
    inputs = n;
    aux = (N-1)/N * ell3*m;     /* Usually we have to send aux, but 1/N times nothing */
    msgs = m + ell3*k;
\\    sizeMPC = inputs + aux + msgs;
    sizeMPC = 972;

    treeM = t*log2(M/t);
    treeN = log2(N);
    bitlen = salt + 
        treeM*seed + /* Initial seeds */
        treeM*hash + /* commitments, Merkle tree */
        t*( treeN*seed + sizeMPC ); 

\\    print("seeds + coms: ", treeM*seed + treeM*hash);
\\    print("t = ", t);
\\    print("per party seeds:", treeN*seed);

    KB = ((bitlen/8)/1000);
    
    return(KB);
}

get_size_mod2mod3_L5(M, N, t)= 
{
    /* Note: t stands for tau, in the write-up's notation */

    local(bitlen, KB, n, k, m, seed, hash, salt, aux, msgs, inputs, treeN, treeM );
  
    /* OWF parameters */
    secpar = 256;       /* Security level (classical; quantum is secpar/2) */
    seed = secpar;      /* Size of a random seed */
    hash = 2*secpar;    /* Digest length for commitments */     
    salt = 256;         /* Bitlength of salt, always 256 */

    sizeMPC = 1943;

    treeM = t*log2(M/t);
    treeN = log2(N);
    bitlen = salt + 
        treeM*seed + /* Initial seeds */
        treeM*hash + /* commitments, Merkle tree */
        t*( treeN*seed + sizeMPC ); 

\\    print("seeds + coms: ", treeM*seed + treeM*hash);
\\    print("t = ", t);
\\    print("per party seeds:", treeN*seed);

    KB = ((bitlen/8)/1000);
    
    return(KB);
}

table6()=
{
    print1("(16, 150, 51) & "); printf("%.02f\n", get_size_mod2mod3(150, 16, 51));
    print1("(16, 168, 45) & "); printf("%.02f\n", get_size_mod2mod3(168, 16, 45));
    print1("(16, 250, 36) & "); printf("%.02f\n", get_size_mod2mod3(250, 16, 36));
    print1("(16, 303, 34) & "); printf("%.02f\n", get_size_mod2mod3(303, 16, 34));
    print("------");
    print1("(64, 151, 45) & "); printf("%.02f\n", get_size_mod2mod3(151, 64, 45));
    print1("(64, 209, 34) & "); printf("%.02f\n", get_size_mod2mod3(209, 64, 34));
    print1("(64, 343, 27) & "); printf("%.02f\n", get_size_mod2mod3(343, 64, 27));
    print1("(64, 515, 24) & "); printf("%.02f\n", get_size_mod2mod3(515, 64, 24));
    print("------");
    print1("(16, 324, 92) & "); printf("%.02f\n", get_size_mod2mod3_L5(324, 16, 92));
    print1("(16, 400, 79) & "); printf("%.02f\n", get_size_mod2mod3_L5(400, 16, 79));
    print1("(16, 604, 68) & "); printf("%.02f\n", get_size_mod2mod3_L5(604, 16, 68));
    print("------");
    print1("(64, 322, 82) & "); printf("%.02f\n", get_size_mod2mod3_L5(322, 64, 82));
    print1("(64, 518, 60) & "); printf("%.02f\n", get_size_mod2mod3_L5(518, 64, 60));
    print1("(64, 604, 57) & "); printf("%.02f\n", get_size_mod2mod3_L5(604, 64, 57));

}

/* max_t = 100 is a good starting point for L1 paramters */
find_params_mod2mod3_L1(max_t, N)=
{
    best_size = 10^20;
    best_M = 0;
    best_t = 0;
    max_M = 3000;
    smallest_t = 10^20;
    
    for(M = 50, max_M, 
        for(t = 8, max_t, 
            if(M-t < 1, 
              break();
             );
            
            if(epsilon(M, N, t) > 128, 

               this_size = get_size_mod2mod3(M, N, t);
               if(t < smallest_t, 
                  print("M = ", M, " t = ", t, " size = ", this_size);
                  smallest_t = t;
               );
               
                if(this_size < best_size, 
                    best_size = this_size;
                    best_M = M;
                    best_t = t;
                );
            );
        );
    );

    print("L1 Done");
    print("best_size = ", best_size);
    print("best_M = ", best_M);
    print("best_t = ", best_t);
    print("N = ", N);
}

/* max_t = 200 is a good starting point for L5 paramters */
find_params_mod2mod3_L5(max_t, N)=
{
    best_size = 10^20;
    best_M = 0;
    best_t = 0;
    max_M = 1000;
    smallest_t = 10^20;
    
    for(M = 50, max_M, 
        for(t = 8, max_t, 
            if(M-t < 1, 
              break();
             );
            
            if(epsilon(M, N, t) > 256, 

               this_size = get_size_mod2mod3_L5(M, N, t);
               if(t < smallest_t, 
                  print("M = ", M, " t = ", t, " size = ", this_size);
                  smallest_t = t;
               );
               
                if(this_size < best_size, 
                    best_size = this_size;
                    best_M = M;
                    best_t = t;
                );
            );
        );
    );

    print("L5 Done");
    print("best_size = ", best_size);
    print("best_M = ", best_M);
    print("best_t = ", best_t);
    print("N = ", N);
}


