function [in] = incorrectValue(in)

while isempty(in) | ~isfloat(in)
    in = input(["Please make in a float, try again for"  ])
end