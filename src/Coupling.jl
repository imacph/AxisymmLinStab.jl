using LegendrePolynomials
using SparseArrays

function Vals_ab_int(vals_a::AbstractVector{Float64},vals_b::AbstractVector{Float64},vals::AbstractVector{Float64},n_x::Int64,x_pts::AbstractVector{Float64})

    int_result::Float64 = 0.

    for k = 1:n_x-1

        int_result += (vals_a[k] * vals[k]* vals_b[k] + vals_a[k+1] * vals[k+1]*vals_b[k+1])*(x_pts[k+1]-x_pts[k])

    end



    return int_result*0.5
end


function calc_coupling_mat_array(n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64,vecs::AbstractVector{Float64})

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)

    u_phi,theta_pts = calc_azim_vel(vecs,n_rad_max,n_leg_max,n_theta_max)
    s_func,theta_pts = calc_stream_func(vecs,n_rad_max,n_leg_max,n_theta_max)
    lsq_func,theta_pts = calc_lsq_stream_func(vecs,n_rad_max,rad_pts,n_leg_max,n_theta_max)



    x_pts = Vector{Float64}(undef,n_theta_max)

    poly_val_mat = Array{Float64,2}(undef,n_theta_max,n_leg_max)
    S_int_val_mat = Array{Float64,2}(undef,n_theta_max,n_leg_max)
    V_int_val_mat = Array{Float64,2}(undef,n_theta_max,n_leg_max)
    T_int_val_mat = Array{Float64,2}(undef,n_theta_max,n_leg_max)

    for k = 1:n_theta_max

        x_pts[k] = cos(theta_pts[n_theta_max-k+1])

        for l = 1:n_leg_max

            poly_val_mat[k,l] = Plm(x_pts[k],l,1)
            V_int_val_mat[k,l] = Plm(x_pts[k],l,1)*x_pts[k]
            T_int_val_mat[k,l] = Plm(x_pts[k],l,1)*(1-x_pts[k]^2)^(1/2)

            if l >= 2
                S_int_val_mat[k,l] = 0.5*(l*(l+1)*Plm(x_pts[k],l,0)-Plm(x_pts[k],l,2))
            else
                S_int_val_mat[k,l] = 0.5*(l*(l+1)*Plm(x_pts[k],l,0))
            end
        end
    end

    norm_facs = zeros(Float64,n_leg_max)
    for l = 1:n_leg_max
        norm_facs[l] = trapezoidal_int(Plm.(x_pts,l,1).^2,n_theta_max,x_pts)
    end


    
    df1_matsp = df1_mat(n_rad_max,rad_pts)
    df2_matsp = df2_mat(n_rad_max,rad_pts)

    n_leg_a::Int64 = 1
    n_leg_b::Int64 = 1

    
    S_ab = zeros(Float64,n_rad_max)
    V_ab = zeros(Float64,n_rad_max)
    T_ab = zeros(Float64,n_rad_max)

    S1_ab = zeros(Float64,n_rad_max)
    V1_ab = zeros(Float64,n_rad_max)
    T1_ab = zeros(Float64,n_rad_max)

    odd_row_odd_col = (n_rad_max-2)*3
    odd_row_even_col = (n_rad_max-2)*3
    
    even_row_odd_col = (n_rad_max-4)*3
    even_row_even_col = (n_rad_max-4)*5

    odd_entries_per_row = (tot_even_leg+1)*odd_row_odd_col + tot_even_leg*odd_row_even_col
    even_entries_per_row = (tot_even_leg+1)*even_row_odd_col + tot_even_leg*even_row_even_col

    tot_entries = even_entries_per_row*(tot_even_leg)+odd_entries_per_row*(tot_even_leg+1)

    mat_array = Array{Float64,2}(undef,tot_entries+1,3)

    oor_vec = zeros(Float64,n_rad_max)

    for k =1:n_rad_max

        oor_vec[k] = 1/rad_pts[k]

    end

    oor_matsp = spdiagm(oor_vec)

    rad_lap_matsp = df2_matsp + *(2*oor_matsp,df1_matsp)
  
    df3_matsp = *(df1_matsp,df2_matsp)
    
    oor_vec = Vector{Float64}(undef,n_rad_max)

    for k = 1:n_rad_max

        oor_vec[k] = 1/rad_pts[k]

    end
    
    oor_matsp = spdiagm(oor_vec)
    oor2_matsp = *(oor_matsp,oor_matsp)
    rad_lap_matsp = df2_matsp + *(2*oor_matsp,df1_matsp) 

    dr_rad_lap_matsp = df3_matsp + *(2*oor_matsp,df2_matsp) + *(-2*oor2_matsp,df1_matsp)
    lap2_matsp = df1_matsp - 2*oor_matsp


    #println(s_func[j,:])
    #println(norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],s_func[j,:],n_theta_max,x_pts))

    for k = 1:tot_even_leg+1
        # odd rows 
        n_leg_a = 2*k-1
        j0 = (n_leg_a-1)*n_rad_max

        norm_fac = (2*n_leg_a+1)/2/(n_leg_a)/(n_leg_a+1)
        norm_fac = 1/norm_facs[n_leg_a]
        #println("odd rows",k)
        for l = 1:tot_even_leg+1
            # odd columns in odd rows 
            
            n_leg_b = 2*l-1
            l0 = (n_leg_b-1)*n_rad_max

            for j = 2:n_rad_max-1
                
                S_ab[j] = norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],s_func[j,:],n_theta_max,x_pts)
                V_ab[j] = norm_fac*Vals_ab_int(V_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],s_func[j,:],n_theta_max,theta_pts)
                T_ab[j] = -norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_a],poly_val_mat[:,n_leg_b],s_func[j,:],n_theta_max,x_pts)- S_ab[j] + V_ab[j]
                #T_ab[j] = norm_fac*Vals_ab_int(T_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],dt_s_func[j,:],n_theta_max,x_pts)
           
            end
            
            dr_S_ab = *(df1_matsp,S_ab)
            dr_V_ab = *(df1_matsp,V_ab)
            
            for j = 2:n_rad_max-1

                i = 3*(j-2) + (l-1)*odd_row_odd_col + (k-1)*odd_entries_per_row

                mat_array[i+1,1] = j + j0
                mat_array[i+1,2] = j-1 + l0
                mat_array[i+1,3] = (V_ab[j]-T_ab[j])/rad_pts[j] * df1_matsp[j,j-1]
                

                mat_array[i+2,1] = j + j0
                mat_array[i+2,2] = j + l0
                mat_array[i+2,3] = (dr_S_ab[j] - dr_V_ab[j] + (S_ab[j] - T_ab[j])/rad_pts[j])/rad_pts[j]
                


                mat_array[i+3,1] = j + j0
                mat_array[i+3,2] = j+1 + l0
                mat_array[i+3,3] = (V_ab[j]-T_ab[j])/rad_pts[j] * df1_matsp[j,j+1]
                
            end

        end

        for l = 1:tot_even_leg
            # even columns in odd rows 

            n_leg_b = 2*l

            l0 = (n_leg_b-1)*n_rad_max

            for j = 2:n_rad_max-1

                S_ab[j] = norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],u_phi[j,:],n_theta_max,x_pts)
                V_ab[j] = norm_fac*Vals_ab_int(V_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],u_phi[j,:],n_theta_max,theta_pts)
                T_ab[j] = -norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_a],poly_val_mat[:,n_leg_b],u_phi[j,:],n_theta_max,x_pts)- S_ab[j] + V_ab[j]
                #T_ab[j] = norm_fac*Vals_ab_int(T_int_val_mat[:,n_leg_b],poly_valmat[:,n_leg_a],dt_u_phi[j,:],n_theta_max,x_pts)
            
            end

            dr_S_ab = *(df1_matsp,S_ab)
            dr_V_ab = *(df1_matsp,V_ab)

            for j = 2:n_rad_max-1

                i = 3*(j-2) + (l-1)*odd_row_even_col + (k-1)*odd_entries_per_row + (tot_even_leg+1)*odd_row_odd_col

                mat_array[i+1,1] = j + j0
                mat_array[i+1,2] = j-1 + l0
                mat_array[i+1,3] = -(V_ab[j]-T_ab[j])/rad_pts[j] * df1_matsp[j,j-1]
                #mat_array[i+1,3] = 0.

                mat_array[i+2,1] = j + j0
                mat_array[i+2,2] = j + l0
                mat_array[i+2,3] = -(dr_S_ab[j] - dr_V_ab[j] + (S_ab[j] - T_ab[j])/rad_pts[j])/rad_pts[j]
                #mat_array[i+2,3] = 0.

                mat_array[i+3,1] = j + j0
                mat_array[i+3,2] = j+1 + l0
                mat_array[i+3,3] = -(V_ab[j]-T_ab[j])/rad_pts[j] * df1_matsp[j,j+1]
                #mat_array[i+3,3] = 0.
            end

        end

    end

    for k = 1:tot_even_leg
        # even rows 
        n_leg_a = 2*k
        j0 = (n_leg_a-1)*n_rad_max
        #println("even rows",k)
        norm_fac = (2*n_leg_a+1)/2/(n_leg_a)/(n_leg_a+1)
        norm_fac = 1/norm_facs[n_leg_a]
        for l = 1:tot_even_leg+1
            # odd columns in even rows 

            n_leg_b = 2*l-1
            l0 = (n_leg_b-1)*n_rad_max
            for j = 2:n_rad_max-1

                S_ab[j] = norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],u_phi[j,:],n_theta_max,x_pts)
                V_ab[j] = norm_fac*Vals_ab_int(V_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],u_phi[j,:],n_theta_max,theta_pts)
                T_ab[j] = -norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_a],poly_val_mat[:,n_leg_b],u_phi[j,:],n_theta_max,x_pts)- S_ab[j] + V_ab[j]
                #T_ab[j] = norm_fac*Vals_ab_int(T_int_val_mat[:,n_leg_b],poly_valmat[:,n_leg_a],dt_u_phi[j,:],n_theta_max,x_pts)
            
            end

            dr_V_ab = *(df1_matsp,V_ab)


            for j = 3:n_rad_max-2

                i = 3*(j-3) + (l-1) * even_row_odd_col + (k-1)*even_entries_per_row + (tot_even_leg+1)*odd_entries_per_row
                
                mat_array[i+1,1] = j+j0
                mat_array[i+1,2] = j-1 +l0
                mat_array[i+1,3] = 2*V_ab[j]/rad_pts[j] * df1_matsp[j,j-1]
                #mat_array[i+1,3] = 0.

                mat_array[i+2,1] = j+j0
                mat_array[i+2,2] = j +l0
                mat_array[i+2,3] = 2*(dr_V_ab[j] + (T_ab[j]+S_ab[j])/rad_pts[j])/rad_pts[j]
                #mat_array[i+2,3]=0.

                mat_array[i+3,1] = j+j0
                mat_array[i+3,2] = j+1 +l0
                mat_array[i+3,3] = 2*V_ab[j]/rad_pts[j] * df1_matsp[j,j+1]
                #mat_array[i+3,3]=0.

            end
        end

        for l = 1:tot_even_leg
            # even columns in even rows 

            n_leg_b = 2*l
            l0 = (n_leg_b-1)*n_rad_max
            for j = 2:n_rad_max-1

                S_ab[j] = norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],s_func[j,:],n_theta_max,x_pts)
                V_ab[j] = norm_fac*Vals_ab_int(V_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],s_func[j,:],n_theta_max,theta_pts)
                T_ab[j] = -norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_a],poly_val_mat[:,n_leg_b],s_func[j,:],n_theta_max,x_pts)- S_ab[j] + V_ab[j]
                #T_ab[j] = norm_fac*Vals_ab_int(T_int_val_mat[:,n_leg_b],poly_valmat[:,n_leg_a],dt_s_func[j,:],n_theta_max,x_pts)

                S1_ab[j] = norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],lsq_func[j,:],n_theta_max,x_pts)
                V1_ab[j] = norm_fac*Vals_ab_int(V_int_val_mat[:,n_leg_b],poly_val_mat[:,n_leg_a],lsq_func[j,:],n_theta_max,theta_pts)
                T1_ab[j] = -norm_fac*Vals_ab_int(S_int_val_mat[:,n_leg_a],poly_val_mat[:,n_leg_b],lsq_func[j,:],n_theta_max,x_pts)- S1_ab[j] + V1_ab[j]
                #T1_ab[j] = norm_fac*Vals_ab_int(T_int_val_mat[:,n_leg_b],poly_valmat[:,n_leg_a],dt_lsq_func[j,:],n_theta_max,x_pts)

            end

            dr_S_ab = *(df1_matsp,S_ab)
            dr_V_ab = *(df1_matsp,V_ab)


            dr_S1_ab = *(df1_matsp,S1_ab)
            dr_V1_ab = *(df1_matsp,V1_ab)

            A = (-dr_S1_ab.+dr_V1_ab)./rad_pts.+(S1_ab.+T1_ab)./rad_pts.^2
            B = (V1_ab.+T1_ab)./rad_pts
            C = (dr_S_ab .+ dr_V_ab)./rad_pts.+(S_ab.+T_ab)./rad_pts.^2
            D = (V_ab.-T_ab)./rad_pts



            for j = 3:n_rad_max-2

                i = 5*(j-3) + (l-1) * even_row_even_col + (k-1)*even_entries_per_row + (tot_even_leg+1)*(odd_entries_per_row + even_row_odd_col)
                
                mat_array[i+1,1] = j+j0
                mat_array[i+1,2] = j-2 +l0
                mat_array[i+1,3] = D[j] * dr_rad_lap_matsp[j,j-2]

                mat_array[i+2,1] = j+j0
                mat_array[i+2,2] = j-1 +l0
                mat_array[i+2,3] = D[j] * (dr_rad_lap_matsp[j,j-1] - n_leg_b*(n_leg_b+1)/rad_pts[j]^2*lap2_matsp[j,j-1])+ B[j] *df1_matsp[j,j-1]+C[j]*rad_lap_matsp[j,j-1]

                mat_array[i+3,1] = j+j0
                mat_array[i+3,2] = j +l0
                mat_array[i+3,3] =  D[j] * (dr_rad_lap_matsp[j,j] - n_leg_b*(n_leg_b+1)/rad_pts[j]^2*lap2_matsp[j,j]) +A[j] +C[j]*(rad_lap_matsp[j,j]-n_leg_b*(n_leg_b+1)/rad_pts[j]^2)

                mat_array[i+4,1] = j+j0
                mat_array[i+4,2] = j+1 +l0
                mat_array[i+4,3] = D[j] * (dr_rad_lap_matsp[j,j+1] - n_leg_b*(n_leg_b+1)/rad_pts[j]^2*lap2_matsp[j,j+1]) + B[j] *df1_matsp[j,j+1]+C[j]*rad_lap_matsp[j,j+1]

                mat_array[i+5,1] = j+j0
                mat_array[i+5,2] = j+2 +l0
                mat_array[i+5,3] =  D[j] * dr_rad_lap_matsp[j,j+2]



            end


        end
        
    end
    mat_array[tot_entries+1,1] = n_rad_max*n_leg_max
    mat_array[tot_entries+1,2] = n_rad_max*n_leg_max
    mat_array[tot_entries+1,3] = 0.
    return mat_array
    
end


function coupling_matrix(n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64,vecs::AbstractVector{ComplexF64})

    mat_array = calc_coupling_mat_array(n_rad_max,rad_ratio,n_leg_max,n_theta_max,real(vecs))

    mat_real = mat_array_to_spmat(mat_array)

    mat_array = calc_coupling_mat_array(n_rad_max,rad_ratio,n_leg_max,n_theta_max,imag(vecs))

    mat_imag = mat_array_to_spmat(mat_array)

    mat = mat_real + 1im*mat_imag

    return mat

end
function save_time_average_azim_vel(n_rad_max::Int64,rad_ratio::Float64,n_leg_max::Int64,n_theta_max::Int64,Ekman::Float64,freq::Float64)

    vecs,mat_G = order_0_parametrized_system(n_rad_max,n_leg_max,rad_ratio,Ekman,freq)

    mat_array = calc_coupling_mat_array(n_rad_max,rad_ratio,n_leg_max,n_theta_max,real(vecs))

    mat_real = mat_array_to_spmat(mat_array)

    mat_array = calc_coupling_mat_array(n_rad_max,rad_ratio,n_leg_max,n_theta_max,imag(vecs))

    mat_imag = mat_array_to_spmat(mat_array)

    mat = mat_real + 1im*mat_imag

    vec_0,mat_G0 = order_0_parametrized_system(n_rad_max,n_leg_max,rad_ratio,Ekman,0.)

    rhs = -0.5*mat*conj.(vecs)

    vecs_0 = mat_G0\rhs

    save_azim_vel(real(vecs_0),n_rad_max,rad_ratio,n_leg_max,n_theta_max)
end


