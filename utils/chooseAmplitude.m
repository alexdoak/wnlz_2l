function in = chooseAmplitude(prompt, amp_max)

    while true

        % Display prompt
        response = input(sprintf( ...
            "%s [enter 'max' for tabletop solitary wave]\n", prompt), 's');

        % Remove whitespace and convert to lowercase
        response = lower(strtrim(response));

        % Maximum amplitude
        if strcmp(response, 'max')
            in = "max";
            break
        end

        % Convert numerical input
        value = str2double(response);

        % Check that input is a valid number
        if isnan(value) || ~isfinite(value)
            fprintf("Please enter a numerical value or 'max'.\n");
            continue
        end

        % Check amplitude range
        if amp_max > 0
            if value < 0 || value > amp_max - 0.01
                fprintf("Amplitude must be between 0 and %g.\n", ...
                    amp_max - 0.01);
                continue
            end

        else
            if value > 0 || value < amp_max + 0.01
                fprintf("Amplitude must be between %g and 0.\n", ...
                    amp_max + 0.01);
                continue
            end
        end

        % Valid amplitude
        in = value;
        break

    end
end