function in = chooseBooleanValue(prompt, default)

    % ---- Default value
    if nargin < 2
        default = [];
    end

    while true

        % ---- Display default in prompt if one exists
        if isempty(default)
            response = input(sprintf("%s\n", prompt), 's');
        else
            if default
                defaultText = "y";
            else
                defaultText = "n";
            end

            response = input(sprintf("%s [default: %s]\n", prompt, defaultText), 's');

            % ---- Empty input -> use default
            if isempty(strtrim(response))
                in = default;
                break
            end
        end

        % ---- Clean up input
        response = lower(strtrim(response));

        % ---- Interpret response
        if any(strcmp(response, {'true', 't', 'yes', 'y', '1'}))
            in = true;
            break

        elseif any(strcmp(response, {'false', 'f', 'no', 'n', '0'}))
            in = false;
            break

        else
            fprintf('Please enter true/false, yes/no, or y/n.\n');
        end

    end
end