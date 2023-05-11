function sdot = gauss_var(t, s)

% [s] = [s(1), s(2), s(3), s(4), s(5), s(6)] 
%     = [  h ,   e ,  f  ,omega,   i ,  w  ]

% parameters
mu = 398600.5;
J2 = 0.00108263;
R = 6378 ;           % Equatorial radius (R/r <1)

h = s(1);  omega = s(4);
e = s(2);      i = s(5);
f = s(3);      w = s(6);

r = h^2/(mu*(1 + e*cos(f)));

p(1) = -(3/2)*(J2*mu*R^2/r^4) * (1 - 3*(sin(i))^2*(sin(w+f))^2);
p(2) = -(3/2)*(J2*mu*R^2/r^4) * (sin(i))^2 * (sin(2*(w+f)));
p(3) = -(3/2)*(J2*mu*R^2/r^4) *  sin(2*i) * sin(w+f);

sdot = [r * p(2);
       (h/mu)*sin(f)*p(1) + (1/(mu*h))*( (h^2 + mu*r)*cos(f) + mu*e*r )*p(2);
       (h/r^2) + (1/(e*h))*( (h^2/mu)*cos(f)*p(1) - (r + h^2/mu)*sin(f)*p(2));
       (r/( h*sin(i) )) * sin(w+f)*p(3)
       (r/h) * cos(w+f)*p(3);
       -(1/(e*h))*( (h^2/mu)*cos(f)*p(1) - (r + h^2/mu)*sin(f)*p(2)) - r*p(3)*(sin(w+f)/(h*tan(i)))];
   

end



