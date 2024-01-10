function A = createCompleteGraph(n)

A = ones(n,n);
for i = 1:n
    A(i,i) = 0;
end

end