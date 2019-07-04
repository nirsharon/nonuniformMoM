function y = myLegendre(l,m,x)
yMat = legendre(l,x);
if m<0
    y = (-1)^abs(m) * factorial(l-abs(m))/factorial(l+abs(m)) * yMat(abs(m)+1,:).';
else
    y = yMat(abs(m)+1,:).';
end

end

