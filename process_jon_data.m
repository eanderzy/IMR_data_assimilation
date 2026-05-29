function bubbles = process_jon_data(file_path,time_centered)
    data = readtable(file_path, 'HeaderLines', 0, 'ReadVariableNames', false);
    data.Properties.VariableNames = {'time', 'radius'};

    bubbles = {};
    current_bubble = [];
    if time_centered == 1
    for i = 1:height(data)
        if data.time(i) == -1000 && data.radius(i) == -1000
            if ~isempty(current_bubble)
                bubbles{end+1} = current_bubble; 
                current_bubble = [];
            end
        else
            current_bubble = [current_bubble; data(i, :)];
        end
    end

    if ~isempty(current_bubble)
        bubbles{end+1} = current_bubble; 
    end
    else
        
        for i = 1:height(data)
        if data.time(i) == -1000 && data.radius(i) == -1000
            if ~isempty(current_bubble)
                min_time = min(current_bubble.time);
                current_bubble.time = current_bubble.time - min_time +5;
                bubbles{end+1} = current_bubble; 
                current_bubble = [];
            end
        else
            current_bubble = [current_bubble; data(i, :)];
        end

        end

    if ~isempty(current_bubble)
        bubbles{end+1} = current_bubble; 
    end
end