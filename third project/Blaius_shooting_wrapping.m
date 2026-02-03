function [actual_ddf_at_zero, counter] = Blaius_shooting_wrapping(~)
    
    %initial_guess = 0.332057336215581;
    initial_guess = rand();
    [actual_ddf_at_zero, counter] = Secant_rf(initial_guess);
end