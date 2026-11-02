function Fhyd = HydroForce(x_dot, y_dot)

    rho = 1025;
    C   = 674;

    vnorm = sqrt(x_dot.^2 + y_dot.^2);    
    scale = -0.5 * rho * C .* vnorm;        

    Fhyd_x = scale .* x_dot;
    Fhyd_y = scale .* y_dot;

    Fhyd = [Fhyd_x(:), Fhyd_y(:)];          % (nT x 2)

end
