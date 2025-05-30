using LegendrePolynomials
using SparseArrays


function calc_theta_grid(n_theta_max::Int64)

    theta_pts = Vector{Float64}(undef,n_theta_max)

    for k = 1:n_theta_max

        theta_pts[k] = k * pi/(n_theta_max+1)
    end

    return theta_pts

end

function calc_cart_grid(n_rad_max::Int64,rad_pts::AbstractVector{Float64},n_theta_max::Int64,theta_pts::AbstractVector{Float64})

    ss = zeros(Float64,(n_rad_max,n_theta_max))
    zz = zeros(Float64,(n_rad_max,n_theta_max))

    for i = 1:n_rad_max
        for j = 1:n_theta_max

            ss[i,j] = rad_pts[i] * sin(theta_pts[j])
            zz[i,j] = rad_pts[i] * cos(theta_pts[j])

        end
    end

    return ss,zz

end


function calc_azim_vel(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,n_leg_max::Int64,n_theta_max::Int64)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    azim_vel_arr = zeros(ComplexF64,(n_rad_max,n_theta_max))

    leg_vec = Vector{ComplexF64}(undef,n_theta_max)

    theta_pts = calc_theta_grid(n_theta_max)
    
    for k = 1:tot_even_leg+1

        n_leg = 2*k-1
        i0 = (2*k-2)*n_rad_max
        
        for j = 1:n_theta_max
            leg_vec[j] = (1+0im)*Plm(cos(theta_pts[j]),n_leg,1)
        end

        azim_vel_arr .+= vecs[i0+1:i0+n_rad_max] * reshape(leg_vec,(1,n_theta_max))
        
    end
    return azim_vel_arr,theta_pts
end

function calc_azim_vel(vecs::AbstractVector{Float64},n_rad_max::Int64,n_leg_max::Int64,n_theta_max::Int64)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    azim_vel_arr = zeros(Float64,(n_rad_max,n_theta_max))

    leg_vec = Vector{Float64}(undef,n_theta_max)

    theta_pts = calc_theta_grid(n_theta_max)
    
    for k = 1:tot_even_leg+1

        n_leg = 2*k-1
        i0 = (2*k-2)*n_rad_max
        
        for j = 1:n_theta_max
            leg_vec[j] = Plm(cos(theta_pts[j]),n_leg,1)
        end

        azim_vel_arr .+= vecs[i0+1:i0+n_rad_max] * reshape(leg_vec,(1,n_theta_max))
        
    end
    return azim_vel_arr,theta_pts
end

function calc_stream_func(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,n_leg_max::Int64,n_theta_max::Int64)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    stream_func_arr = zeros(ComplexF64,(n_rad_max,n_theta_max))

    leg_vec = Vector{ComplexF64}(undef,n_theta_max)

    theta_pts = calc_theta_grid(n_theta_max)
    
    for k = 1:tot_even_leg

        n_leg = 2*k
        i0 = (2*k-1)*n_rad_max
        
        for j = 1:n_theta_max
            leg_vec[j] = (1+0im)*Plm(cos(theta_pts[n_theta_max-j+1]),n_leg,1)
        end

        stream_func_arr .+= vecs[i0+1:i0+n_rad_max] * reshape(leg_vec,(1,n_theta_max))
        
    end
    return stream_func_arr,theta_pts
end

function calc_stream_func(vecs::AbstractVector{Float64},n_rad_max::Int64,n_leg_max::Int64,n_theta_max::Int64)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    stream_func_arr = zeros(Float64,(n_rad_max,n_theta_max))

    leg_vec = Vector{Float64}(undef,n_theta_max)

    theta_pts = calc_theta_grid(n_theta_max)
    
    for k = 1:tot_even_leg

        n_leg = 2*k
        i0 = (2*k-1)*n_rad_max
        
        for j = 1:n_theta_max
            leg_vec[j] = Plm(cos(theta_pts[n_theta_max-j+1]),n_leg,1)
        end

        stream_func_arr .+= vecs[i0+1:i0+n_rad_max] * reshape(leg_vec,(1,n_theta_max))
        
    end
    return stream_func_arr,theta_pts
end

function calc_lsq_stream_func(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,rad_pts::AbstractVector{Float64},n_leg_max::Int64,n_theta_max::Int64)

    df1_matsp = df1_mat(n_rad_max,rad_pts)
    df2_matsp = df2_mat(n_rad_max,rad_pts)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    lsq_stream_func_arr = zeros(ComplexF64,(n_rad_max,n_theta_max))

    leg_vec = Vector{ComplexF64}(undef,n_theta_max)

    theta_pts = calc_theta_grid(n_theta_max)

    oor_vec = Vector{Float64}(undef,n_rad_max)

    for j = 1:n_rad_max
        oor_vec[j] = 1/rad_pts[j]
    end

    oor_matsp = spdiagm(oor_vec)
    oor2_matsp = *(oor_matsp,oor_matsp)


    rad_lap_matsp = df2_matsp + *(2*oor_matsp,df1_matsp)

    for k = 1:tot_even_leg

        n_leg = 2*k
        i0 = (2*k-1)*n_rad_max
        
        for j = 1:n_theta_max
            leg_vec[j] = (1+0im)*Plm(cos(theta_pts[n_theta_max-j+1]),n_leg,1)
        end

        lsq_stream_func_arr .+=  *(rad_lap_matsp-n_leg*(n_leg+1)*oor2_matsp,vecs[i0+1:i0+n_rad_max]) * reshape(leg_vec,(1,n_theta_max))
        
    end
    return lsq_stream_func_arr,theta_pts
end


function calc_kinetic_energy(vecs::AbstractVector{Float64},n_rad_max::Int64,n_leg_max::Int64,n_theta_max::Int64,rad_ratio::Float64)

    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    
    v_phi,theta_pts = calc_azim_vel(vecs,n_rad_max,n_leg_max,n_theta_max)
    v_theta,theta_pts = calc_colat_vel(vecs,n_rad_max,rad_pts,n_leg_max,n_theta_max)
    v_r,theta_pts = calc_rad_vel(vecs,n_rad_max,rad_pts,n_leg_max,n_theta_max)
    

    kin_energy = 0.5*(v_phi.^2+v_theta.^2+v_r.^2)

    kin_energy_radii = zeros(Float64,n_rad_max)

    for j =1:n_rad_max

        kin_energy_radii[j] = trapezoidal_int(kin_energy[j,:].*sin.(theta_pts),n_theta_max,theta_pts)
    end

    kin_energy_tot = trapezoidal_int(kin_energy_radii.*rad_pts.^2,n_rad_max,rad_pts)

    return kin_energy_tot
    
end



function calc_colat_vel(vecs::AbstractVector{Float64},n_rad_max::Int64,rad_pts::AbstractVector{Float64},n_leg_max::Int64,n_theta_max::Int64)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    df1_matsp = df1_mat(n_rad_max,rad_pts)

    colat_vel_arr = zeros(Float64,(n_rad_max,n_theta_max))

    leg_vec = Vector{Float64}(undef,n_theta_max)

    theta_pts = calc_theta_grid(n_theta_max)
    
    for k = 1:tot_even_leg

        n_leg = 2*k
        i0 = (2*k-1)*n_rad_max
        
        for j = 1:n_theta_max
            leg_vec[j] = Plm(cos(theta_pts[n_theta_max-j+1]),n_leg,1)
        end

        colat_vel_arr .+= (*(df1_matsp,vecs[i0+1:i0+n_rad_max])-vecs[i0+1:i0+n_rad_max]./rad_pts) * reshape(leg_vec,(1,n_theta_max))
        
    end
    return colat_vel_arr,theta_pts
end

function calc_rad_vel(vecs::AbstractVector{Float64},n_rad_max::Int64,rad_pts::AbstractVector{Float64},n_leg_max::Int64,n_theta_max::Int64)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    
    rad_vel_arr = zeros(Float64,(n_rad_max,n_theta_max))

    leg_vec = Vector{Float64}(undef,n_theta_max)


    theta_pts = calc_theta_grid(n_theta_max)
    
    for k = 1:tot_even_leg

        n_leg = 2*k
        i0 = (2*k-1)*n_rad_max
        
        for j = 1:n_theta_max
            x = cos(theta_pts[n_theta_max-j+1])
            leg_vec[j] = Plm(x,n_leg,1)*x/(1-x^2)^(1/2)-x*dnPl(x,n_leg,1)+(1-x^2)*dnPl(x,n_leg,2)
        end

        rad_vel_arr .+= (vecs[i0+1:i0+n_rad_max]./rad_pts) * reshape(leg_vec,(1,n_theta_max))

        
    end
    return rad_vel_arr,theta_pts
end

function calc_lsq_stream_func(vecs::AbstractVector{Float64},n_rad_max::Int64,rad_pts::AbstractVector{Float64},n_leg_max::Int64,n_theta_max::Int64)

    df1_matsp = df1_mat(n_rad_max,rad_pts)
    df2_matsp = df2_mat(n_rad_max,rad_pts)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    lsq_stream_func_arr = zeros(Float64,(n_rad_max,n_theta_max))

    leg_vec = Vector{Float64}(undef,n_theta_max)

    theta_pts = calc_theta_grid(n_theta_max)

    oor_vec = Vector{Float64}(undef,n_rad_max)

    for j = 1:n_rad_max
        oor_vec[j] = 1/rad_pts[j]
    end

    oor_matsp = spdiagm(oor_vec)
    oor2_matsp = *(oor_matsp,oor_matsp)


    rad_lap_matsp = df2_matsp + *(2*oor_matsp,df1_matsp)

    for k = 1:tot_even_leg

        n_leg = 2*k
        i0 = (2*k-1)*n_rad_max
        
        for j = 1:n_theta_max
            leg_vec[j] = Plm(cos(theta_pts[n_theta_max-j+1]),n_leg,1)
        end

        lsq_stream_func_arr .+=  *(rad_lap_matsp-n_leg*(n_leg+1)*oor2_matsp,vecs[i0+1:i0+n_rad_max]) * reshape(leg_vec,(1,n_theta_max))
        
    end
    return lsq_stream_func_arr,theta_pts
end

function save_azim_vel(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64)

    azim_vel_arr,theta_pts = calc_azim_vel(vecs,n_rad_max,n_leg_max,n_theta_max)
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm("v_phi_real.txt",real(azim_vel_arr))
    writedlm("v_phi_imag.txt",imag(azim_vel_arr))

    writedlm("ss.txt",ss)
    writedlm("zz.txt",zz)

    nothing
end

function save_azim_vel(vecs::AbstractVector{Float64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64)

    azim_vel_arr,theta_pts = calc_azim_vel(vecs,n_rad_max,n_leg_max,n_theta_max)
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm("v_phi_real.txt",real(azim_vel_arr))
    writedlm("v_phi_imag.txt",imag(azim_vel_arr))

    writedlm("ss.txt",ss)
    writedlm("zz.txt",zz)

    nothing
end

function save_lsq_stream_func(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64)
    
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    lsq_stream_func_arr,theta_pts = calc_lsq_stream_func(vecs,n_rad_max,rad_pts,n_leg_max,n_theta_max)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm("lsq_func_real.txt",real(lsq_stream_func_arr))
    writedlm("lsq_func_imag.txt",imag(lsq_stream_func_arr))

    writedlm("ss.txt",ss)
    writedlm("zz.txt",zz)

end

function save_stream_func(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64)

    stream_func_arr,theta_pts = calc_stream_func(vecs,n_rad_max,n_leg_max,n_theta_max)
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm("s_func_real.txt",real(stream_func_arr))
    writedlm("s_func_imag.txt",imag(stream_func_arr))

    writedlm("ss.txt",ss)
    writedlm("zz.txt",zz)

    nothing
end



function save_stream_func(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64, directory::String)

    
    stream_func_arr,theta_pts = calc_stream_func(vecs,n_rad_max,n_leg_max,n_theta_max)
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm(directory*"/s_func_real.txt",real(stream_func_arr))
    writedlm(directory*"/s_func_imag.txt",imag(stream_func_arr))

    writedlm(directory*"/ss.txt",ss)
    writedlm(directory*"/zz.txt",zz)

    nothing
end

function save_stream_func(vecs::AbstractVector{Float64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64, directory::String)

    
    stream_func_arr,theta_pts = calc_stream_func(vecs,n_rad_max,n_leg_max,n_theta_max)
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm(directory*"/s_func_real.txt",real(stream_func_arr))
    writedlm(directory*"/s_func_imag.txt",imag(stream_func_arr))

    writedlm(directory*"/ss.txt",ss)
    writedlm(directory*"/zz.txt",zz)

    nothing
end

function save_azim_vel(vecs::AbstractVector{ComplexF64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64, directory::String)

    
    azim_vel_arr,theta_pts = calc_azim_vel(vecs,n_rad_max,n_leg_max,n_theta_max)
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm(directory*"/v_phi_real.txt",real(azim_vel_arr))
    writedlm(directory*"/v_phi_imag.txt",imag(azim_vel_arr))

    writedlm(directory*"/ss.txt",ss)
    writedlm(directory*"/zz.txt",zz)

    nothing
end

function save_azim_vel(vecs::AbstractVector{Float64},n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64, directory::String)
    

    azim_vel_arr,theta_pts = calc_azim_vel(vecs,n_rad_max,n_leg_max,n_theta_max)
    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)
    ss,zz = calc_cart_grid(n_rad_max,rad_pts,n_theta_max,theta_pts)

    writedlm(directory*"/v_phi_real.txt",real(azim_vel_arr))
    writedlm(directory*"/v_phi_imag.txt",imag(azim_vel_arr))

    writedlm(directory*"/ss.txt",ss)
    writedlm(directory*"/zz.txt",zz)

    nothing
end