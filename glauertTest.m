% Right now, this is just a test.  This basically copies
% the problem from the 12 week exam as a proof of concept.
function[] = glauertTest()
    b = 10;
    function[x] = a0(y)
        x = 2*pi*ones(1, numel(y));
    end
    function[x] = alpha0(y)
        x = zeros(1, numel(y));
    end
    function[x] = c(y)
        x = ones(1, numel(y));
    end
    A = glauertAn(0.1, @a0, @alpha0, @c, b, [0 2.5]);
    % A should be [0.01536 0.00140] to match the values
    % for the 12 week exam (values for AOA == 0.1 rad)
    display(A);
end