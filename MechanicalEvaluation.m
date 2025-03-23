 
filename='Table.xlsx';

% read from table
A=readtable(filename,'sheet','Tabelle1');
file=table2array(A(:,1));
area=table2array(A(:,2));
einsp=table2array(A(:,3));
ticker=table2array(A(:,4));

% Initialize arrays for calculations
maxforce=zeros(length(file),1);
maxstress=zeros(length(file),1);
maxstrain_displ=zeros(length(file),1);
maxstrain=zeros(length(file),1);
totaldispl=zeros(length(file),1);
Emodul=zeros(length(file),1);


% Placeholder for grouping
sth='nix';

% Plot handles and Labels for Dropdown-Menu
plotHandles = [];
plotLabels = {};

% Creating a figure for the plots
fig = figure();

% Loop for processing each file
for i=1:length(file)
    
    if strcmp(ticker(i),sth)==0
        if i > 1
            createDropdownMenu(plotHandles, plotLabels);
            plotHandles = [];
            plotLabels = {};
        end
        figure();
        hold on;
        sth=ticker(i);
    end
    
    % Read raw data from the txt. files
    C=readtable(string(file(i)));
    % Smoothing the force values according to the Gauss function
    force=table2array(smoothdata(C(:,2),'gaussian',50)); 
    displacement=table2array(C(:,3)); 
    
    % calculation for Stress and Strain
    stress=force/area(i);
    strain=displacement/einsp(i);




    
    % calculation of max. Stress 
    maxstress(i)=max(stress);
    maxforce(i)=max(force);
    for a=1:length(stress)
        if stress(a)==maxstress(i)
            break
        end
    end

    % Cut off curve after maximum
    if length(stress)<a
        stressCut=stress;
        strainCut=strain;
    else
        stressCut=stress(1:a);
        strainCut=strain(1:a);
    end


    % calculation fof max. Strain
    maxstrain(i) = max(strainCut);


    %%%%%% Cut off lead time and set curves to zero point %%%%%

    % Lead time correction: Define threshold value (Very small value to eliminate initial noise)
    threshold = 0.00001; 
    startIndex = find(stressCut >= threshold, 1);

    % If a valid starting point was found, cut the lead time
    if ~isempty(startIndex)
        stressCut = stressCut(startIndex:end);
        strainCut = strainCut(startIndex:end);
    end

    % Zero point correction for all curves, all curves start at the origin
    if ~isempty(stressCut) && ~isempty(strainCut)
        stressCut = stressCut - stressCut(1); 
        strainCut = strainCut - strainCut(1); 
    end


    %%%%%% Calculate total displacement after adjustment %%%%%%

    % Deformation margin after correction, from origin (0/0) to break    
    if ~isempty(strainCut)
        totaldispl(i) = strainCut(end) - strainCut(1);  
    end    

    

    %%%%%% Calculate Young's modulus %%%%%%
    
    % Parameter for the window size (number of points for the local regression) 
    window_size = 20;  % Can be adjusted (smaller values = more sensitive to noise)
    
    % Calculation of local gradients with sliding window
    num_points = length(strainCut);
    slopes = zeros(num_points - window_size, 1);
    
    for j = 1:(num_points - window_size)
        % Select the current window for the linear regression
        strain_window = strainCut(j:j+window_size);
        stress_window = stressCut(j:j+window_size);
        
        % Linear regression for the window
        p = polyfit(strain_window, stress_window, 1);
        
        % Save the calculated gradient
        slopes(j) = p(1);
    end
    
    % Find the highest and most stable gradient (youngs-modulus)
    [max_slope, max_idx] = max(slopes);
    
    % Take the mean value over a small window around the maximum gradient
    stable_range = max(max_idx-5,1):min(max_idx+5,length(slopes));
    Emodul(i) = mean(slopes(stable_range)); 




    %%%%%% Plotting the curves %%%%%%

    % Plotting the adjusted stress-strain curve and creating the drop-down menu
    h = plot(strainCut, stressCut*1000, 'DisplayName', file{i},'LineWidth', 1);
    plotHandles = [plotHandles, h];
    plotLabels{end+1} = file{i};
    
    % Set axis labeling & title
    xlabel('strain [-]', 'FontSize', 12);
    ylabel('stress [GPa]', 'FontSize', 12);
    title(sth, 'FontSize', 10);
    set(gca, 'FontSize', 14);  

end


%%%%%% Fill table with calculated values %%%%%%

T = table(maxforce, maxstress*1000, totaldispl, maxstrain,Emodul*1000,'VariableNames', {'MaxForce [N]','MaxStress [GPa]','TotalDispl [mm]','MaxStrain [-]','Emodul [GPa]'});
writetable(T,filename,'sheet','Tabelle1','Range','E:I');



%%%%%% create the boxplot %%%%%%

A1 = 1000 * Emodul;  
names = "1.Version";    % Initialization of the names for the groups (ticker categories)
anz = length(names);    % Number of groups for the boxplot
xax = 1;                % Control of the X-axis labeling (1 = show, 0 = hide)

ticker{end+1} = 'last';
sth = ticker(1);
start = 1;

figure()                
b = 0;          % group counter        


    % Define the color sheme
    colors = [
        179/255,     242/255,     253/255
        132/255,     239/255,     236/255
         99/255,     231/255,     207/255
         98/255,     220/255,     169/255
        112/255,     209/255,     131/255
        128/255,     197/255,      95/255
        143/255,     182/255,      60/255
        154/255,     161/255,      34/255
        157/255,     138/255,      28/255
        156/255,     117/255,      36/255
        153/255,      99/255,      48/255
        151/255,      82/255,      59/255
        149/255,      65/255,      71/255
        147/255,      48/255,      84/255
        144/255,      30/255,      99/255
        140/255,       2/255,     115/255

    ];

    % Run through all lines of ticker (grouping variable)
    for i = 1:length(ticker(:,1))

        % Check whether a new group starts (i.e. 'ticker' changes)
        if strcmp(ticker(i), sth) == 0
            sth = ticker(i);
            ende = i - 1;

            % Update the group counter
            b = b + 1;
            names{b} = ticker{i-1};
            B = A1(start:ende, 1:anz);

            % Create X-values for the boxplot (group on the X-axis)
            X = ones(ende-start+1, 1) * b;

            hold on

            
            %%% Create boxplot manually %%%

            % Prepare data
            data = B(:,1);
            q1 = quantile(data, 0.25); % first quartile
            q3 = quantile(data, 0.75); % third quartile
            iqr = q3 - q1;             % Interquartile distance
            median_val = median(data); % Median
            lower_whisker = max(min(data), q1 - 1.5 * iqr);
            upper_whisker = min(max(data), q3 + 1.5 * iqr);

            % draw box
            hold on;
            fill([b-0.3, b+0.3, b+0.3, b-0.3], [q1, q1, q3, q3], colors(mod(b-1, size(colors, 1)) + 1, :), 'EdgeColor', 'k', 'LineWidth', 2);

            % draw median line
            line([b-0.3, b+0.3], [median_val, median_val], 'Color', 'k', 'LineWidth', 2);

            % draw whiskers 
            line([b, b], [lower_whisker, q1], 'Color', 'k', 'LineWidth', 2);
            line([b, b], [q3, upper_whisker], 'Color', 'k', 'LineWidth', 2);

            % Cross lines at the whisker ends
            line([b-0.1, b+0.1], [lower_whisker, lower_whisker], 'Color', 'k', 'LineWidth', 2); 
            line([b-0.1, b+0.1], [upper_whisker, upper_whisker], 'Color', 'k', 'LineWidth', 2); 

            % draw Scatter-Points         
            scatter(b * ones(size(data)), data, 40, 'k', 'filled');

            
            % Set the new start index for the next group
            start = i;
        end
    end

    % Axis labeling and formatting
    if xax == 1
        set(gca, "XTick", [1:b], "XTickLabel", names, "fontsize", 24); 
        xtickangle(45);  % Oblique angle of the X-axis labels
    else
        set(gca, "XTick", []);
    end

    ylabel('E-Modul [GPa]', 'FontSize', 24);  
    hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%% Create dropdown menu %%%%%

function createDropdownMenu(handles, labels)
    if isempty(handles)
        return;
    end
    
    fig = gcf;
    menuItems = strcat(labels, ' (Visible)');
    menu = uicontrol('Style', 'popupmenu', ...
                     'String', menuItems, ...
                     'Units', 'normalized', ...
                     'Position', [0.8 0 0.15 0.05], ...  % Position outside the plot area
                     'Callback', @(src, event) toggleVisibility(src, handles, labels));
end

% Callback function to toggle visibility
function toggleVisibility(src, handles, labels)
    val = src.Value;
    if strcmp(handles(val).Visible, 'on')
        handles(val).Visible = 'off';
    else
        handles(val).Visible = 'on';
    end
    updateMenuLabels(src, handles, labels);
end

% Update menu labels to reflect visibility status
function updateMenuLabels(src, handles, labels)
    newString = labels;
    for i = 1:length(handles)
        if strcmp(handles(i).Visible, 'on')
            newString{i} = strcat(labels{i}, ' (Visible)');
        else
            newString{i} = strcat(labels{i}, ' (Hidden)');
        end
    end
    src.String = newString;
end



