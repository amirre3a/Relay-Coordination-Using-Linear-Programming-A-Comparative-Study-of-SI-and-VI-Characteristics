% EX3 Amirreza Shafieinejad $402449105
clc;
clear;

% 1. Constants
n_relays = 14;        % Number of relays
CTI = 0.2;            % Critical Time Interval

% Pickup currents (Ip)
Ip = [600; 800; 500; 800; 600; 500; 600; 500; 600; 500; 600; 500; 600; 800];

% Fault currents (If_main, If_backup)
If_main = [3232; 5924; 5924; 3556; 3783; 2401; 6109; 6093; 2484; 3883; 3707; 5899; 2991; 5199];
If_backup = [3232; 996; 1890; 3556; 3783; 2401; 1874; 1890; 2484; 2344; 3707; 987; 2991; 996];
If_backup = max(If_backup, Ip * 1.1); % Ensure If_backup >= 1.1 * Ip

% Standard Inverse
a_standard = 0.14;
b_standard = 0.02;
% Very Inverse
a_very = 13.5;
b_very = 1.0;       

% TMS bounds
TMS_min = 0.05;
TMS_max = 2;

% optimization variable TMS
TMS = optimvar('TMS', n_relays, 'LowerBound', TMS_min, 'UpperBound', TMS_max);

% Instructions for constructing A
instructions = [
    1 1 1 -1;   1 6 6  1;
    2 1 1  1;    2 2 2 -1;
    3 2 2 -1;    3 7 7  1;
    4 2 2  1;    4 3 3 -1;
    5 3 3  1;    5 4 4 -1;
    6 4 4  1;    6 5 5 -1;
    7 6 6 -1;    7 14 14 1;
    8 7 7 -1;    8 13 13 1;
    9 7 7  1;    9 8 8 -1;
    10 9 9 -1;   10 10 10 1;
    11 10 10 -1; 11 11 11 1;
    12 11 11 -1; 12 12 12 1;
    13 12 12 -1; 13 13 13 1;
    14 12 12 -1; 14 14 14 1;
    15 8 8 1;    15 3 3 -1;
    16 1 1 1;    16 14 14 -1
];

n_constraints = max(instructions(:,1));
B = -CTI * ones(n_constraints, 1);
slack = optimvar('slack', n_constraints, 'LowerBound', 0);

% Solve Standard Inverse
[TMS_SI, slack_SI] = scenarioSolveAndRecord('report.txt', TMS, Ip, If_main, a_standard, b_standard, instructions, B, slack, 'Standard Inverse', If_main, If_backup);

% Solve Very Inverse
[TMS_VI, slack_VI] = scenarioSolveAndRecord('report.txt', TMS, Ip, If_main, a_very, b_very, instructions, B, slack, 'Very Inverse', If_main, If_backup);

% Plot comparison if both scenarios succeeded
if ~isempty(TMS_SI) && ~isempty(TMS_VI)
    figure;
    subplot(2,1,1);
    bar([TMS_SI, TMS_VI]);
    title('TMS Comparison: Standard Inverse vs Very Inverse');
    xlabel('Relay Number'); ylabel('TMS');
    legend('Standard Inverse','Very Inverse'); grid on;

    subplot(2,1,2);
    bar([slack_SI, slack_VI]);
    title('Slack Comparison');
    xlabel('Constraint Number'); ylabel('Slack');
    legend('Standard Inverse Slack','Very Inverse Slack'); grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%  Function  %%%%%%%%%%%%%%%%%%%%%%%
function [TMS_sol, slack_sol] = scenarioSolveAndRecord(reportFile, TMS, Ip, If_main, a, b, instructions, B, slack, scenarioName, If_main_vector, If_backup_vector)


n_relays = length(Ip);
n_constraints = max(instructions(:,1));

% Compute M
M = a ./ ((If_main ./ Ip).^b - 1);

% Build A
A = zeros(n_constraints, n_relays);
for i = 1:size(instructions,1)
    r = instructions(i,1);
    c = instructions(i,2);
    m_idx = instructions(i,3);
    s = instructions(i,4);
    A(r,c) = s * M(m_idx);
end

% Objective and constraints
objective = sum(TMS .* (a ./ ((If_main ./ Ip).^b - 1))) + 1e6 * sum(slack);
constraints = A * TMS - slack <= B;

problem = optimproblem('Objective', objective, 'Constraints', constraints);
[sol, ~, exitflag] = solve(problem);

fprintf('--- %s Results ---\n', scenarioName);
if exitflag == 1
    TMS_sol = sol.TMS;
    slack_sol = sol.slack;
    fprintf('Optimal TMS values:\n');
    disp(TMS_sol);
    fprintf('Slack values:\n');
    disp(slack_sol);
    avg_TMS = mean(TMS_sol);
    fprintf('Average TMS: %.4f\n', avg_TMS);

    main_times = TMS_sol .* (a ./ ((If_main_vector ./ Ip).^b - 1));
    backup_times = TMS_sol .* (a ./ ((If_backup_vector ./ Ip).^b - 1));

    fid = fopen(reportFile, 'a');
    if fid == -1
        error('Cannot open report file.');
    end

    fprintf(fid, '--- %s Run ---\n', scenarioName);
    fprintf(fid, 'Relay\tTMS\tMainTime(s)\tBackupTime(s)\n');
    for rr = 1:n_relays
        fprintf(fid, '%d\t%.4f\t%.4f\t%.4f\n', rr, TMS_sol(rr), main_times(rr), backup_times(rr));
    end
    fprintf(fid, 'Average TMS: %.4f\n\n', avg_TMS);
    fclose(fid);

else
    TMS_sol = []; slack_sol = [];
    fprintf('No feasible solution found for %s. Exit flag: %d\n', scenarioName, exitflag);
end
fprintf('\n');
end
