function [u_diff,v_diff] = diffusion(Nx,gamma,del_x,del_y,u_cell_sides_nrdr,v_cell_sides_nrdr,u_cell_tb_nrdr,v_cell_tb_nrdr)
    u_diff = 0*ones(Nx+2,Nx+2);
    v_diff = 0*ones(Nx+2,Nx+2);
    V = del_x*del_y;
    for i = 2 : Nx + 1
        for j = 2 : Nx + 1
            u_R1 = u_cell_sides_nrdr(i-1,j) * del_y;
            u_R2 = - u_cell_sides_nrdr(i-1,j-1) * del_y;
            u_R3 = u_cell_tb_nrdr(i-1,j-1) * del_x;
            u_R4 = - u_cell_tb_nrdr(i,j-1) * del_x;
            u_diff(i,j) = gamma * (u_R1 + u_R2 + u_R3 + u_R4)/V;
    
            v_R1 = v_cell_sides_nrdr(i-1,j) * del_y;
            v_R2 = - v_cell_sides_nrdr(i-1,j-1) * del_y;
            v_R3 = v_cell_tb_nrdr(i-1,j-1) * del_x;
            v_R4 = - v_cell_tb_nrdr(i,j-1) * del_x;
            v_diff(i,j) = gamma * (v_R1 + v_R2 + v_R3 + v_R4)/V;
        end
    end
end