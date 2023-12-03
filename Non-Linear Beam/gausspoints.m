function [w, gp] = gausspoints(ng)

if ng == 1
    w = 2;
    gp = 0.00;
elseif ng == 2
    w = [1.00 1.00];
    gp = [1/sqrt(3) -1/sqrt(3)];
elseif ng == 3
    w = [5/9 8/9 5/9];
    gp = [-sqrt(3/5) 0.00 sqrt(3/5)];
elseif ng ==4
w = [0.6521451548 0.6521451548 0.3478548451 0.3478548451];     
gp = [0.3399810435 -0.3399810435 0.8611363116 -0.8611363116]; 
end    