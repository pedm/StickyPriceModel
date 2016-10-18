if sum(imag(x_start_unbounded)) ~= 0
    x_start_unbounded
    error('Error: your bounds do not make sense. You are using an initial value that lies outside the bounds.')
end