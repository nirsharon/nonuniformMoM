function Ysh = mySph_v2(l,m,theta,phi)
% Complex valued Spherical Harmonics
% Ysh = sqrt((2*l+1)*factorial(l-m)./4./pi./factorial(l+m)) .* myLegendre(l,m,cos(theta)) .* exp(1i*m*phi);
[theta_unique,~,theta_map] = unique(theta);
[phi_unique,~,phi_map] = unique(phi);
Ysh = bsxfun(@times,sqrt((2*l+1)*factorial(l-m)./4./pi./factorial(l+m)), myLegendre(l,m,cos(theta_unique)));
W = exp(1i*phi_unique*m);
Ysh = Ysh(theta_map,:) .* W(phi_map,:);

% Real valued Spherical Harmonics
% if m>0
%     Ysh = sqrt(2)*sqrt((2*l+1)*factorial(l-m)/4/pi/factorial(l+m)) .* myLegendre(l,m,cos(theta)) .* cos(m*phi);
% elseif m==0
%     Ysh = sqrt((2*l+1)*factorial(l-m)/4/pi/factorial(l+m)) .* myLegendre(l,m,cos(theta));
% else
%     Ysh = sqrt(2)*sqrt((2*l+1)*factorial(l-abs(m))/4/pi/factorial(l+abs(m))) .* myLegendre(l,abs(m),cos(theta)) .* sin(abs(m)*phi);
% end

end

