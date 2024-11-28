function NetworkSimulation(initial_values_csv, n, t)
    %%%Simulation of a M*N neural network, where initial_value is a path to
    %%%a .csv file containing initial parameters for neurons in grid 
    
    
    
    [neuron_state, initial_values] = calculate_display_graph(initial_values_csv, n, 0);
    for i = 1:t
        new_amyloid = move_time_step(neuron_state, initial_values, n, i*12);
        calculate_display_graph(new_amyloid, n, i);
    end 
    
    
    
end     
function ab_t = amyloid_time(t, y, ab_0)
    %%%returns amyloid at time t normalized for grid using y as as the rate
    %%%constant and ab_0 is the initial amyloid con 
    dxdt = @(t, ab) y*ab - y*ab^2;
    tspan = [0 t];
    [T, ab] = ode45(dxdt, tspan, ab_0);
    ab_t = ab(end);
end 
function new_amyloid_values_csv = move_time_step(neuron_state, initial_values, n, t)
    new_values = zeros(n, 1);
    disp(new_values);
    neuron_state = reshape(neuron_state, 1, []);
    for i = 1:n
        ab = amyloid_time(t, 0.08, initial_values(i));
        if (neuron_state(1, i)==1)
            ab = (1-0.7)*ab;
        end 
        if (neuron_state(1, i)==2)
            ab = (1-0.3)*ab;
        end 
        new_values (i) = ab;
    end 
    new_amyloid_values_csv = sprintf('amyloid_time_%d.csv', t);
    disp(new_values);
    writematrix(new_values, new_amyloid_values_csv); 
end
    
function [neuron_state, initial_values] = calculate_display_graph(initial_values_csv, n, t)
    %%%Introduce arrays
    m=sqrt(n);
    neuron_state = zeros(m,m);
    initial_values=readmatrix(initial_values_csv);
    %disp(initial_values);
    
    
    %%%Set initial neuron values 
    blocked_synapses = floor(initial_values/0.125);
    blocked_synapses = reshape(blocked_synapses, [m,m])';
    disp(blocked_synapses);
    
    %%%Calculate neuron_state 
    blocked_synapses_total = blocked_synapses;
    for i = 1:m
        for j = 1:m
            to_add = 0; 
            if (i>1) 
                %blocked_synapses_total(i,j) = blocked_synapses_total(i,j)+blocked_synapses(i-1,j);
                if blocked_synapses(i-1, j)> 4
                    to_add = to_add + mod (blocked_synapses(i-1, j), 4);
                end 
                if (j>1&&blocked_synapses(i-1, j-1)>4)
                    to_add = to_add +1;
                end 
                if (j<m &&blocked_synapses(i-1, j+1)>6)
                    to_add = to_add +1;
                end 
            end 
            if (j>1)
                %blocked_synapses_total(i,j) = blocked_synapses_total(i,j)+blocked_synapses(i,j-1);
                if (blocked_synapses(i, j-1) > 2 && blocked_synapses(i, j-1) < 6)
                    to_add = to_add + blocked_synapses(i, j-1) - 2; 
                end 
            end 
            if (i<m) 
                %blocked_synapses_total(i,j) = blocked_synapses_total(i,j)+blocked_synapses(i+1,j);
                if (blocked_synapses(i+1, j) > 3)
                    to_add = to_add + 3; 
                else 
                    if(blocked_synapses(i+1, j)>0)
                        to_add = to_add + blocked_synapses(i+1, j);
                    end
                end 
                if(j>1&&blocked_synapses(i+1, j-1)>2)
                    to_add = to_add + 1;
                end 
                if(j<m && blocked_synapses(i+1, j+1)>0)
                    to_add = to_add +1;
                end 
            end 
            if (j<m)
                %blocked_synapses_total(i,j) = blocked_synapses_total(i,j)+blocked_synapses(i,j+1);
                if (blocked_synapses(i, j+1) > 1)
                    to_add = to_add + 1; 
                end 
                if (blocked_synapses(i, j+1) > 6)
                    to_add = to_add + mod(blocked_synapses(i, j+1), 6); 
                end
            end
           
            blocked_synapses_total(i,j) = blocked_synapses_total(i,j) + to_add; 
            
            if (i==1&&(j==1||j==m)||i==m&&(j==1||j==m)) %%corner
                if (blocked_synapses_total(i,j)<=2)
                    neuron_state(i,j)=0; 
                else 
                    if (blocked_synapses_total(i,j)<=4)
                        neuron_state(i,j)=1;
                    else 
                        if (blocked_synapses_total(i,j)<8)
                            neuron_state(i,j)= 2;
                        else 
                            neuron_state(i,j) = 3;
                        end 
                    end 
                end
            else  
    
                if (i==1||j==m||i==m||j==1) %%row//column
                    if (blocked_synapses_total(i,j)<=2)
                        neuron_state(i,j)=0; 
                    else 
                        if (blocked_synapses_total(i,j)<=5)
                            neuron_state(i,j)=1; 
                        else 
                            if (blocked_synapses_total(i,j)<=10)
                                neuron_state(i,j)= 2; 
                            else 
                                neuron_state(i,j) = 3; 
                            end 
                        end 
                    end
                else 
                    if (blocked_synapses_total(i,j)<=3)
                        neuron_state(i,j)=0; 
                    else 
                        if (blocked_synapses_total(i,j)<=6)
                            neuron_state(i,j)=1; 
                        else 
                            if (blocked_synapses_total(i,j)<=12)
                                neuron_state(i,j)= 2; 
                            else 
                                neuron_state(i,j) = 3;  
                            end 
                        end 
                    end
                end 
            end 
            if(blocked_synapses(i,j) == 8) %%instant death if all synapses blocked
                neuron_state(i,j) = 3; 
            end
        end 
    end 
    disp(neuron_state)
    grid = zeros(3*m, 3*m+1); 
    for i = 1:m
        for j = 1:m
            row_start = 3*(i-1) + 1;
            col_start = 3*(j-1) + 1;
            grid (row_start +1, col_start+1) = neuron_state(i, j);
            nb_to_block = blocked_synapses(i,j);
            clockwise_order = [
                row_start,     col_start;     % Top-left
                row_start,     col_start+1;   % Top-center
                row_start,     col_start+2;   % Top-right
                row_start+1,   col_start+2;   % Middle-right
                row_start+2,   col_start+2;   % Bottom-right
                row_start+2,   col_start+1;   % Bottom-center
                row_start+2,   col_start;     % Bottom-left
                row_start+1,   col_start;     % Middle-left
            ];
            for k = 1:8
                if k <= nb_to_block
                    % Blocked cell (purple = 5)
                    grid(clockwise_order(k, 1), clockwise_order(k, 2)) = 5;
                else
                    % Unblocked cell (blue = 4)
                    grid(clockwise_order(k, 1), clockwise_order(k, 2)) = 4;
                end
            end
        end
    end
    grid(1, 3*m+1) = 0; 
    grid(2, 3*m+1) = 1;
    grid(3, 3*m+1) = 2;
    grid(4, 3*m+1) = 3;
    grid(5, 3*m+1) = 4;
    grid(6, 3*m+1) = 5;
    imagesc(grid);
    colormap([
        %1, 1, 1; % 0 = white (background)
        1, 0.75, 0.8; % 0 = pink
        0, 1, 0; % 1 = green
        1, 0, 0; % 2 = red
        0, 0, 0; % 3 = black
        0, 0, 1; % 4 = blue (unblocked)
        0.5, 0, 0.5; % 5 = purple (blocked)
    ]);
    %colormap(turbo);
    axis equal tight;
    title('Neurons and Synapses');
    to_save = sprintf('figure%i.jpg', t);
    saveas(gcf, to_save); % Save as MATLAB .fig file
   
        
end
            
    

    
    
    