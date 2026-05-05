function [u_cell_sides_nrdr,v_cell_sides_nrdr,u_cell_tb_nrdr,v_cell_tb_nrdr] = nrml_drv(Nx,u_cell_cntr,v_cell_cntr,del_x,del_y)
    u_cell_sides_nrdr = zeros(Nx,Nx+1);
    v_cell_sides_nrdr = zeros(Nx,Nx+1);
    
    u_cell_tb_nrdr = zeros(Nx+1,Nx);
    v_cell_tb_nrdr = zeros(Nx+1,Nx);
    
    for i = 2 : Nx+1
        for j = 1 : Nx+1
            u_cell_sides_nrdr(i-1,j) = (u_cell_cntr(i,j+1) - u_cell_cntr(i,j))/del_x ; 
            v_cell_sides_nrdr(i-1,j) = (v_cell_cntr(i,j+1) - v_cell_cntr(i,j))/del_x ; 
        end
    end
    
    for k = 1 : Nx+1
        for l = 2 : Nx+1
            u_cell_tb_nrdr(k,l-1) = (u_cell_cntr(k,l) - u_cell_cntr(k+1,l))/del_y;
            v_cell_tb_nrdr(k,l-1) = (v_cell_cntr(k,l) - v_cell_cntr(k+1,l))/del_y;
        end
    end
end