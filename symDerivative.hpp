double dfdp(double FAO, double g, double p, double Del, double pp);
double dfdpp(double p, double Del, double pp);
double dfdg(double p, double FAO, double g);
double dfdFAO(double g, double p);
double d2fdp2(double FAO, double g, double p, double Del, double pp);
double d2fdpp2(double p, double Del, double pp);
double d2fdg2(double p, double FAO, double g);
double psi_prime_dot0(double FAO, double g, double p, double Del, double pp);
double psi_prime_ddot0(double FAO, double g, double p, double Del, double pp);
void sig_generation(double sig1, double sig2,
                    double sig_dot1, double sig_dot2,
                    double sig_ddot1, double sig_ddot2);
double sigProd(double sig[1][2]);
void sigProd_generation(double sig1, double sig2);
double sigProdA_dot(double sigA_dot1, double sigA_dot2, double sigA_dot3,
                    double sigA_dot4, double sigA_dot5, double sigA_dot6,
                    double sigA1, double sigA2, double sigA3,
                    double sigA4, double sigA5, double sigA6);
double sigProdD_dot(double sigD_dot1, double sigD_dot2, double sigD_dot3, double sigD_dot4,
                    double sigD1, double sigD2, double sigD3, double sigD4);
double d2Beta_barF_dBeta2(double beta, double a0, double b0, double n0);