function X = gradient_descent(cost_function, grad_function, initial_guess, weight,plot_flag)
    
    learning_rate = 1;
    cost = abs(cost_function(initial_guess,weight));
    last_cost = cost + 1000000000000000;
    Cost = [];
    while(abs(abs(cost) - abs(last_cost))/abs(last_cost) >= 0.0001 && learning_rate >= 0.00000001)
        X = initial_guess - learning_rate * grad_function(initial_guess,weight);
        step_cost = abs(cost_function(X,weight));
        if step_cost >= cost
            learning_rate = learning_rate * 0.5;
        else
            last_cost = cost;
            cost = step_cost;
            initial_guess = X;
            learning_rate = learning_rate * 1.1;
            Cost = [Cost , cost];
        end
    end
    if plot_flag
        plot(1:size(Cost,2),Cost);
        hold off;
    end
end