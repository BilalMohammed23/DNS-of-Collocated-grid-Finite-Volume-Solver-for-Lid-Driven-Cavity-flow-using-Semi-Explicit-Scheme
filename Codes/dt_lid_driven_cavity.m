function [delta_t_max] = dt_lid_driven_cavity(u_cell_cntr,v_cell_cntr,del_x,del_y,Nx,gamma)
    CFL = 0.5;
    
    umax = max(abs(u_cell_cntr(2:Nx+1,2:Nx+1)),[],'all');
    vmax = max(abs(v_cell_cntr(2:Nx+1,2:Nx+1)),[],'all');
    
    % avoid division by zero at initial condition
    if umax < 1e-12
        umax = 1e-12;
    end
    
    if vmax < 1e-12
        vmax = 1e-12;
    end
    
    dt_adv = 1 / (umax/del_x + vmax/del_y);
    
    dt_diff = 1 / (2*gamma*(1/del_x^2 + 1/del_y^2));
    
    delta_t_max = CFL * min(dt_adv, dt_diff);
end