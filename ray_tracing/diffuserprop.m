function xout = diffuserprop(xin, thetain, Dx, z0, n, np)

    xout = xin + z0 * (n*thetain/np + (n/np-1)*Dx);
    
