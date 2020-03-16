

    

    if (norm == "W")
    {
        double scal_V = 0.0;
        double scal_H = 0.0;

        


        for (int i = 1; i < M; i++)//1...N
            scal_H += (Ro_0(h * i, 1) - H[i]) * (Ro_0(h * i, 1) - H[i]);
            
        for (int i = 1; i < M; i++)//1...N
            scal_V += (u_0 (h * i, 1) - V[i]) * (u_0 (h * i, 1) - V[i]);
        
        scal_H *= h;
        scal_V *= h;
        
        scal_H += 0.5 * h * ((Ro_0(0, 1) - H[0]) * (Ro_0(0, 1) - H[0])
                             + (Ro_0(h * M, 1) - H[M]) * (Ro_0(h * M, 1) - H[M])
                            );
        scal_V += 0.5 * h * ((u_0 (h * 0, 1) - V[0]) * (u_0 (h * 0, 1) - V[0])
                             + (u_0 (h * M, 1) - V[M]) * (u_0 (h * M, 1) - V[M]));

        resid.first = scal_H;
        resid.second = scal_V;
        
        double temp1 = 0.;
        for (int i = 0; i < M; i++)
        {
            temp1 += (u_0 (h * (i), 1) - V[i] - u_0 (h * (i + 1), 1) + V[i + 1]) / h * 
                     (u_0 (h * (i), 1) - V[i] - u_0 (h * (i + 1), 1) + V[i + 1]) / h;
        }
        temp1 *= h;

        double temp2 = 0.;
        for (int i = 0; i < M; i++)
        {
            temp2 += (Ro_0 (h * (i), 1) - H[i] + H[i + 1] - Ro_0 (h * (i + 1), 1)) / h * 
                     (Ro_0 (h * (i), 1) - H[i] + H[i + 1] - Ro_0 (h * (i + 1), 1)) / h;
        }
        temp2 *= h;
        resid.first = sqrt(fabs (resid.first + temp2));
        resid.second = sqrt(fabs (resid.second + temp1));
        


    }
