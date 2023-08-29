function y = mmka(S, Km, n, eta)
 % positive regulation
 
y = (S^n + eta) / (S^n + Km^n);

end 