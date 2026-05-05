function [u_conv,v_conv] = convection(Nx,u_cell_cntr,v_cell_cntr,u_cell_face_sides_intp,v_cell_face_tb_intp,del_x,del_y)
    u_conv = 0*ones(Nx+2,Nx+2);
    v_conv = 0*ones(Nx+2,Nx+2);
    V = del_x*del_y;
    for i = 2 : Nx+1
        for j = 2 : Nx+1
            u_R1 = ((u_cell_cntr(i,j) + u_cell_cntr(i,j+1))/2) * u_cell_face_sides_intp(i-1,j+1) *del_y; 
            u_R2 = ((u_cell_cntr(i,j) + u_cell_cntr(i,j-1))/2) * u_cell_face_sides_intp(i-1,j) * del_y;
            u_R3 = ((u_cell_cntr(i,j) + u_cell_cntr(i-1,j))/2) * v_cell_face_tb_intp(i,j-1) * del_x;
            u_R4 = ((u_cell_cntr(i,j) + u_cell_cntr(i+1,j))/2) * v_cell_face_tb_intp(i+1,j-1) * del_x;
            u_conv(i,j) = (u_R1 - u_R2 + u_R3 - u_R4)/V;
    
            v_R1 = ((v_cell_cntr(i,j) + v_cell_cntr(i,j+1))/2) * u_cell_face_sides_intp(i-1,j+1) *del_y; 
            v_R2 = ((v_cell_cntr(i,j) + v_cell_cntr(i,j-1))/2) * u_cell_face_sides_intp(i-1,j) * del_y;
            v_R3 = ((v_cell_cntr(i,j) + v_cell_cntr(i-1,j))/2) * v_cell_face_tb_intp(i,j-1) * del_x;
            v_R4 = ((v_cell_cntr(i,j) + v_cell_cntr(i+1,j))/2) * v_cell_face_tb_intp(i+1,j-1) * del_x;
            v_conv(i,j) = (v_R1 - v_R2 + v_R3 - v_R4)/V;
        end
    end
end