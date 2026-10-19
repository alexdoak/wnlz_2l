function in = chooseWaveDirection(prompt, default)

    % Default value
    if nargin < 2
        default = [];
    end

    while true

        % Display prompt
        if isempty(default)
            response = input(sprintf("%s\n", prompt), 's');
        else
            response = input(sprintf("%s [default: %s]\n", prompt, default), 's');

            % Empty input -> use default
            if isempty(strtrim(response))
                in = default;
                break
            end
        end

        % Clean up input
        response = lower(strtrim(response));

        % Interpret response
        if any(strcmp(response, {'r', 'right'}))
            in = "right";
            break

        elseif any(strcmp(response, {'l', 'left'}))
            in = "left";
            break

        else
            fprintf("Please enter 'r' or 'l'.\n");
        end

    end
end