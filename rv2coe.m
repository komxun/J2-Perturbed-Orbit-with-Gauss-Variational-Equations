function [a, e_norm, i, omega, w, f, h_norm] = rv2coe(r,v,mu)

h = cross(r,v);    % Angular momentum vector
h_norm = norm(h);
k = [0; 0; 1];     % Unit vector pointing to the normal vector of the plane
n = cross(k, h);   % Nodal vector

e = (-cross(h,v)/mu) -(r/norm(r));
e_norm = norm(e);                        % eccentricity    
a = (h_norm^2)/(mu*(1-e_norm^2));      % semi-major axis
i = acos(dot(k,h)/norm(h));              % inclination

%% Right ascension of the ascending node, Omega
if (n(2) < 0)
    omega = 2*pi - acos(n(1)/norm(n));
else
    omega = acos(n(1)/norm(n));
end
omega_deg = omega * 180/pi;

%% Argument of periaposis, w
if (dot(e,k)<0)
    w = 2*pi - acos(dot(n,e)/(norm(n)*norm(e)));
else
    w = acos(dot(n,e)/(norm(n)*norm(e)));
end
w_deg = w*180/pi;

%% True anomaly, f
if (dot(r,v)<0)
    f = 2*pi - acos(dot(r,e)/(norm(r)*norm(e)));
else
    f = acos(dot(r,e)/(norm(r)*norm(e)));
end
f_deg = f*180/pi;

end


