function D = discrete_derivative_segment(n, h)
    d = [-1 1; -1 1];
    D = 1/h * kron(speye(n), d);
end