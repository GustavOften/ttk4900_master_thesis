function [delta,d_delta,dd_delta,ddd_delta] = delta_5_cosine(phi)
        phi = mod(phi,2*pi);
        a =      0.1148;
        b =    -0.04023;
        c =   -0.001535;
        d =    0.001152;
        e =   0.0008215;
        delta =  a+b*cos(2*phi)+c*cos(4*phi)+d*cos(6*phi)+e*cos(8*phi);
        d_delta =  -b*2*sin(2*phi)-c*4*sin(4*phi)-d*6*sin(6*phi)-e*8*sin(8*phi);
        dd_delta =  -b*4*cos(2*phi)-c*16*cos(4*phi)-d*36*cos(6*phi)-e*64*cos(8*phi);
        ddd_delta =  b*8*sin(2*phi)+c*16*4*sin(4*phi)+d*36*6*sin(6*phi)+e*64*8*sin(8*phi);
end

