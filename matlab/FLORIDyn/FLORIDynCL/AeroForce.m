function Faero = AeroForce(FT, yaw_iT)
    Fx =  FT .* cos(yaw_iT);
    Fy = -FT .* sin(yaw_iT);
    Faero = [Fx, Fy];

end
