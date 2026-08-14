function r = order_parameter( ang, NN )

NN = size(ang,2);
r = (1/NN) * sum( exp(1i*ang), 2 );
r = abs(r);
