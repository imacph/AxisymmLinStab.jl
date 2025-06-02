

function prep_time_step(n_rad_max::Int64,n_leg_max::Int64,rad_ratio::Float64,Ekman::Float64,freq::Float64)

    if iseven(n_leg_max)
        n_leg_max += 1
    end

    diffusion_matsp,coriolis_matsp,time_matsp,bc_matsp = spatial_mats(n_rad_max,n_leg_max,rad_ratio)

    mat_G1 = 1im*freq*time_matsp + (1+0im)*(coriolis_matsp-Ekman*diffusion_matsp+bc_matsp)
    mat_G0 = coriolis_matsp-Ekman*diffusion_matsp+bc_matsp
    mat_G2 = 2im*freq*time_matsp + (1+0im)*(coriolis_matsp-Ekman*diffusion_matsp+bc_matsp)

    rhs1 = zeros(ComplexF64,n_rad_max*n_leg_max)
    rhs1[n_rad_max] = -1.

    vecs_1 = mat_G1\rhs1

    coupling_mat1 = coupling_matrix(n_rad_max,rad_ratio,n_leg_max,n_leg_max,vecs_1)

    rhs0 = -0.5*coupling_mat1 * conj(vecs_1)
    rhs2 = -0.5*coupling_mat1 * vecs_1

    vecs_0 = mat_G0\rhs0
    vecs_2 = mat_G2\rhs2
    coupling_mat2 = coupling_matrix(n_rad_max,rad_ratio,n_leg_max,n_leg_max,vecs_2)
    coupling_mat0 = coupling_matrix(n_rad_max,rad_ratio,n_leg_max,n_leg_max,vecs_0)

    return diffusion_matsp,coriolis_matsp,time_matsp,bc_matsp,coupling_mat1,coupling_mat0,coupling_mat2,Ekman,freq,n_rad_max,n_leg_max,rad_ratio

end


function time_step(dt::Float64,n_dt::Int64,Re::Float64,prep_input::Tuple)

    diffusion_matsp,coriolis_matsp,time_matsp,bc_matsp,coupling_mat1,coupling_mat0,coupling_mat2,Ekman,freq,n_rad_max,n_leg_max,rad_ratio = prep_input

    vec_t0 = zeros(Float64,n_rad_max*n_leg_max)
    vec_t1 = zeros(Float64,n_rad_max*n_leg_max)
    vec_t2 = zeros(Float64,n_rad_max*n_leg_max)

    t::Float64 = 0

    epsi::Float64 = Re * Ekman^(1/2)
    
    mat_G = 1.5*time_matsp + dt*(coriolis_matsp-Ekman*diffusion_matsp) + bc_matsp

    vec_ave = zeros(Float64,n_rad_max*n_leg_max)
    vec_t0[n_rad_max-1] = epsi
    vec_t1[n_rad_max-1] = epsi

    KE_array = zeros(Float64,n_dt)

    precond = factorize(mat_G)
    for i=1:n_dt
        
        mat_val = real(exp(1im*freq*t)*coupling_mat1+exp(2im*freq*t)*coupling_mat2*epsi+coupling_mat0*epsi)*epsi

        #vec_t2 = ( mat_G+dt*mat_val) \ (2*time_matsp*vec_t1-0.5*time_matsp*vec_t0)
        vec_t2 = bicgstabl!(vec_t2,mat_G+dt*mat_val,2*time_matsp*vec_t1-0.5*time_matsp*vec_t0,Pl=precond)

        
        vec_t0 = vec_t1
        
        vec_t1 = vec_t2
        
        t += dt
        vec_ave += vec_t2/n_dt
        
        
        KE_array[i] = calc_kinetic_energy(vec_t2, n_rad_max, n_leg_max, n_leg_max, rad_ratio)
        
        
    end


    
    return vec_ave,KE_array
end


function time_step(min_dt::Float64,n_dt::Int64,Re::Float64,prep_input::Tuple,to::TimerOutput)

    dt = min_dt
    diffusion_matsp,coriolis_matsp,time_matsp,bc_matsp,coupling_mat1,coupling_mat0,coupling_mat2,Ekman,freq,n_rad_max,n_leg_max,rad_ratio = prep_input

    vec_t0 = zeros(Float64,n_rad_max*n_leg_max)
    vec_t1 = zeros(Float64,n_rad_max*n_leg_max)
    vec_t2 = zeros(Float64,n_rad_max*n_leg_max)

    t::Float64 = 0

    epsi::Float64 = Re * Ekman^(1/2)
    
    mat_G = 1.5*time_matsp + dt*(coriolis_matsp-Ekman*diffusion_matsp) + bc_matsp

    vec_ave = zeros(Float64,n_rad_max*n_leg_max)
    vec_t0[n_rad_max-1] = epsi
    vec_t1[n_rad_max-1] = epsi

    KE_array = zeros(Float64,n_dt)

    precond = lu(mat_G)

    coupling_mat1_re = real(coupling_mat1)*epsi
    coupling_mat1_im = imag(coupling_mat1)*epsi
    
    coupling_mat2_re = real(coupling_mat2)*epsi^2
    coupling_mat2_im = imag(coupling_mat2)*epsi^2

    coupling_mat0_re = real(coupling_mat0)*epsi^2
    coupling_mat0_im = imag(coupling_mat0)*epsi^2
    
    # (A+iB) * (C+iD) = (AC-BD) + i(AD+BC)

    for i=1:n_dt
        
        @timeit to "coupling update" (mat_val = (cos(freq*t)*coupling_mat1_re .- sin(freq*t)*coupling_mat1_im .+
                    cos(2*freq*t)*coupling_mat2_re .- sin(2*freq*t)*coupling_mat2_im .+
                    coupling_mat0_re .- coupling_mat0_im ))
        #   mat_val = real(exp(1im*freq*t)*coupling_mat1+exp(2im*freq*t)*coupling_mat2*epsi+coupling_mat0*epsi)*epsi
        

        #@timeit to "time step" (vec_t2 = ( mat_G+dt*mat_val) \ (2*time_matsp*vec_t1-0.5*time_matsp*vec_t0))
        @timeit to "time step" (vec_t2 = bicgstabl!(vec_t2,mat_G+dt*mat_val,2*time_matsp*vec_t1-0.5*time_matsp*vec_t0,Pl=precond))

        
        vec_t0 = vec_t1
        
        vec_t1 = vec_t2
        
        t += dt
        vec_ave += vec_t2/n_dt
        
        
        @timeit to "calc KE" (KE_array[i] = calc_kinetic_energy(vec_t2, n_rad_max, n_leg_max, n_leg_max, rad_ratio))
        
        #dt = min(sqrt(2*KE_array[i])/(4/3*pi*(1-rad_ratio^3)/(1-rad_ratio)^3))
        #println("Time step: ", dt)
        #println((sqrt(2*KE_array[i])/(4/3*pi*(1-rad_ratio^3)/(1-rad_ratio)^3)))
    end


    
    return vec_ave,KE_array
end

