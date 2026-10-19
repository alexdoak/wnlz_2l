function in = chooseIntegerValue(prompt, minval, maxval, default)

    % ---- Default limits
    if nargin < 2 || isempty(minval)
        minval = -inf;
    end

    if nargin < 3 || isempty(maxval)
        maxval = inf;
    end

    % ---- Change order of minval-maxval if the bounds are wrong way round
    if minval > maxval
        temp = minval;
        minval = maxval;
        maxval = temp;
    end

    % ---- Default value
    if nargin < 4
        default = [];
    end

    while true

        if isempty(default)
            in = input(sprintf("%s\n", prompt));
        else
            in = input(sprintf("%s [default value: %d]\n", prompt, default));

            % ---- Empty input -> use default
            if isempty(in)
                in = default;
            end
        end

        % ---- Check that input is a scalar integer
        if ~isnumeric(in) || ~isscalar(in) || ~isfinite(in) || in ~= fix(in)
            fprintf('Please enter a single integer value.\n');
            continue
        end

        % ---- Check minimum
        if in <= minval
            fprintf('Must be greater than %g, try again.\n', minval);
            continue
        end

        % ---- Check maximum
        if in >= maxval
            fprintf('Must be less than %g, try again.\n', maxval);
            continue
        end

        % ---- Passed all checks
        break

    end
end