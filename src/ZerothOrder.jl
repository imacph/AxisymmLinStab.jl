using SparseArrays
using DelimitedFiles

function calc_rad_pts(n_rad_max::Int64,rad_ratio::Float64)

    r_o::Float64 = 1/(1-rad_ratio)
    r_i::Float64 = r_o * rad_ratio

    rad_pts = Vector{Float64}(undef,n_rad_max)
    for k = 1:n_rad_max
        rad_pts[k] = 0.5*(cos((n_rad_max-k)*pi/(n_rad_max-1))+r_o+r_i)
    end

    return r_i,r_o,rad_pts
end


function loop_df1_array!(df1_array::AbstractArray{Float64,2},r_1::Float64,r_2::Float64,k::Int64)

    df1_array[2*k+1,1] = k+1
    df1_array[2*k+1,2] = k
    df1_array[2*k+1,3] = -1/(r_2-r_1)

    df1_array[2*k+2,1] = k+1
    df1_array[2*k+2,2] = k+2
    df1_array[2*k+2,3] = -df1_array[2*k+1,3]
    
    nothing
end


function loop_df2_array!(df2_array::AbstractArray{Float64,2},n_rad_max::Int64,rad_pts::AbstractVector{Float64})

    df2_array[1,1] = 1
    df2_array[1,2] = 1
    df2_array[1,3] = -2/(rad_pts[2]-rad_pts[1])*(rad_pts[3]-rad_pts[1])
    
    df2_array[2,1] = 1
    df2_array[2,2] = 2
    df2_array[2,3] = 2/(rad_pts[2]-rad_pts[1])/(rad_pts[2]-rad_pts[3])

    df2_array[3,1] = 1
    df2_array[3,2] = 3
    df2_array[3,3] = 2/(rad_pts[3]-rad_pts[1])/(rad_pts[2]-rad_pts[3])

    df2_array[3*n_rad_max-2,1] = n_rad_max
    df2_array[3*n_rad_max-2,2] = n_rad_max
    df2_array[3*n_rad_max-2,3] = -2/(rad_pts[n_rad_max-1]-rad_pts[n_rad_max])*(rad_pts[n_rad_max-2]-rad_pts[n_rad_max])

    df2_array[3*n_rad_max-1,1] = n_rad_max
    df2_array[3*n_rad_max-1,2] = n_rad_max-1
    df2_array[3*n_rad_max-1,3] = 2/(rad_pts[n_rad_max-1]-rad_pts[n_rad_max])/(rad_pts[n_rad_max-1]-rad_pts[n_rad_max-2])

    df2_array[3*n_rad_max,1] = n_rad_max
    df2_array[3*n_rad_max,2] = n_rad_max-2
    df2_array[3*n_rad_max,3] = 2/(rad_pts[n_rad_max-2]-rad_pts[n_rad_max])/(rad_pts[n_rad_max-1]-rad_pts[n_rad_max-2])

    for k = 1:(n_rad_max-2)

        df2_array[3*k+1,1] = k+1
        df2_array[3*k+1,2] = k
        df2_array[3*k+1,3] = 2/(rad_pts[k+1]-rad_pts[k])/(rad_pts[k+2]-rad_pts[k])

        df2_array[3*k+2,1] = k+1
        df2_array[3*k+2,2] = k+1
        df2_array[3*k+2,3] = -2/(rad_pts[k+2]-rad_pts[k+1])/(rad_pts[k+1]-rad_pts[k])

        df2_array[3*k+3,1] = k+1
        df2_array[3*k+3,2] = k+2
        df2_array[3*k+3,3] = 2/(rad_pts[k+2]-rad_pts[k+1])/(rad_pts[k+2]-rad_pts[k])
    end
    nothing
end

function df1_mat(n_rad_max::Int64,rad_pts::AbstractVector{Float64})


    df1_array = Array{Float64,2}(undef,2*n_rad_max,3)

    df1_array[1,1] = 1
    df1_array[1,2] = 1
    df1_array[1,3] = -1/(rad_pts[2] - rad_pts[1])

    df1_array[2,1] = 1
    df1_array[2,2] = 2
    df1_array[2,3] = -df1_array[1,3]

    df1_array[2*n_rad_max-1,1] = n_rad_max
    df1_array[2*n_rad_max-1,2] = n_rad_max-1
    df1_array[2*n_rad_max-1,3] = -1/(rad_pts[n_rad_max] - rad_pts[n_rad_max-1])

    df1_array[2*n_rad_max,1] = n_rad_max
    df1_array[2*n_rad_max,2] = n_rad_max
    df1_array[2*n_rad_max,3] = -df1_array[2*n_rad_max-1,3]

    for k = 1:(n_rad_max-2)

        r_2 = rad_pts[k+2]
        r_1 = rad_pts[k]
        loop_df1_array!(df1_array,r_1,r_2,k)
    end

    df1_spmat = sparse(df1_array[:,1],df1_array[:,2],df1_array[:,3])

    return df1_spmat
end

function df2_mat(n_rad_max::Int64,rad_pts::AbstractVector{Float64})

    
    df2_array = Array{Float64,2}(undef,3*n_rad_max,3)

    loop_df2_array!(df2_array,n_rad_max,rad_pts)

    
    df2_spmat = sparse(df2_array[:,1],df2_array[:,2],df2_array[:,3])

    return df2_spmat
end

function calc_diff_mat_array(n_rad_max::Int64,tot_even_leg::Int64,rad_pts::AbstractVector{Float64},df1_matsp::AbstractArray{Float64},df2_matsp::AbstractArray{Float64},df3_matsp::AbstractArray{Float64},df4_matsp::AbstractArray{Float64})

    tot_even_entries = ((n_rad_max-4)*5+4)*tot_even_leg
    tot_odd_entries = ((n_rad_max-2)*3+2)*(tot_even_leg+1) 
    
    even_row_tot = (n_rad_max-4)*5+4
    odd_row_tot = (n_rad_max-2)*3+2
    comb_row_tot = even_row_tot + odd_row_tot

    mat_array = Array{Float64,2}(undef,tot_even_entries+tot_odd_entries,3)



    mat_array[1,1] = 1
    mat_array[1,2] = 1
    mat_array[1,3] = 0.

    mat_array[2,1] = n_rad_max
    mat_array[2,2] = n_rad_max
    mat_array[2,3] = 0.

    for j = 2:n_rad_max-1

        i = 3*(j-1)

        mat_array[i,1] = j
        mat_array[i,2] = j-1
        mat_array[i,3] = df2_matsp[j,j-1] + 2/rad_pts[j] * df1_matsp[j,j-1]

        mat_array[i+1,1] = j
        mat_array[i+1,2] = j
        mat_array[i+1,3] = df2_matsp[j,j] - 2/(rad_pts[j]*rad_pts[j])
        
        mat_array[i+2,1] = j
        mat_array[i+2,2] = j+1
        mat_array[i+2,3] = df2_matsp[j,j+1] + 2/rad_pts[j] * df1_matsp[j,j+1]

    end
        

    for k = 1:tot_even_leg

        n_leg_even = 2*k
        n_leg_odd = 2*k+1
        
        i = (k-1)*comb_row_tot+odd_row_tot

        j0 = (2*k-1)*n_rad_max

        mat_array[i+1,1] = j0+1
        mat_array[i+1,2] = j0+1
        mat_array[i+1,3] = 0.

        mat_array[i+2,1] = j0+2
        mat_array[i+2,2] = j0+2
        mat_array[i+2,3] = 0.

        mat_array[i+3,1] = j0 + n_rad_max-1
        mat_array[i+3,2] = j0 + n_rad_max-1
        mat_array[i+3,3] = 0.
        
        mat_array[i+4,1] = j0 + n_rad_max
        mat_array[i+4,2] = j0 + n_rad_max
        mat_array[i+4,3] = 0.

        for j = 3:n_rad_max-2

            i = 5*(j-2)+(k-1)*comb_row_tot+odd_row_tot
            
            mat_array[i,1] = j + j0
            mat_array[i,2] = j-2 + j0
            mat_array[i,3] = df4_matsp[j,j-2] +4/rad_pts[j] * df3_matsp[j,j-2]
            
            mat_array[i+1,1] = j+ j0
            mat_array[i+1,2] = j-1+ j0
            mat_array[i+1,3] = df4_matsp[j,j-1] +4/rad_pts[j] * df3_matsp[j,j-1] - 2*n_leg_even*(n_leg_even+1)/(rad_pts[j]*rad_pts[j])*df2_matsp[j,j-1] 

            mat_array[i+2,1] = j+ j0
            mat_array[i+2,2] = j+ j0
            mat_array[i+2,3] = df4_matsp[j,j] +4/rad_pts[j] * df3_matsp[j,j]- 2*n_leg_even*(n_leg_even+1)/(rad_pts[j]*rad_pts[j])*df2_matsp[j,j] + n_leg_even*(n_leg_even^3+2*n_leg_even^2-n_leg_even-2)/rad_pts[j]^4

            mat_array[i+3,1] = j+ j0
            mat_array[i+3,2] = j+1+ j0
            mat_array[i+3,3] = df4_matsp[j,j+1] +4/rad_pts[j] * df3_matsp[j,j+1] - 2*n_leg_even*(n_leg_even+1)/(rad_pts[j]*rad_pts[j])*df2_matsp[j,j+1] 

            mat_array[i+4,1] = j+ j0
            mat_array[i+4,2] = j+2+ j0
            mat_array[i+4,3] = df4_matsp[j,j+2] +4/rad_pts[j] * df3_matsp[j,j+2]

        end

        i = k*comb_row_tot

        j0 = 2*k*n_rad_max
        mat_array[i+1,1] = j0+1
        mat_array[i+1,2] = j0+1
        mat_array[i+1,3] = 1.

        mat_array[i+2,1] = j0+n_rad_max
        mat_array[i+2,2] = j0+n_rad_max
        mat_array[i+2,3] = 1.

        for j = 2:n_rad_max-1

            i = 3*(j-1)+k*comb_row_tot

            mat_array[i,1] = j+j0
            mat_array[i,2] = j-1+j0
            mat_array[i,3] = df2_matsp[j,j-1] + 2/rad_pts[j] * df1_matsp[j,j-1]

            mat_array[i+1,1] = j+j0
            mat_array[i+1,2] = j+j0
            mat_array[i+1,3] = df2_matsp[j,j] - n_leg_odd*(n_leg_odd+1)/(rad_pts[j]*rad_pts[j])
            
            mat_array[i+2,1] = j+j0
            mat_array[i+2,2] = j+1+j0
            mat_array[i+2,3] = df2_matsp[j,j+1] + 2/rad_pts[j] * df1_matsp[j,j+1]

        end


    end

    return mat_array
end

function calc_cor_mat_array(n_rad_max::Int64,tot_even_leg::Int64,rad_pts::AbstractVector{Float64},df1_matsp::AbstractArray{Float64})


    even_entries_per_block = (n_rad_max-4)*3
    even_entries_per_row::Int64 = 2*even_entries_per_block

    odd_entries_per_block::Int64 = (n_rad_max-2)*3
    odd_entries_per_row::Int64 = 2*odd_entries_per_block

    mat_array = Array{Float64,2}(undef,even_entries_per_row*tot_even_leg+odd_entries_per_row*tot_even_leg+1,3)
    
    for j = 2:n_rad_max-1

        i = 3*(j-2)

        mat_array[i+1,1] = j
        mat_array[i+1,2] = j-1 + n_rad_max
        mat_array[i+1,3] = -2 *3/5*df1_matsp[j,j-1]

        mat_array[i+2,1] = j
        mat_array[i+2,2] = j + n_rad_max
        mat_array[i+2,3] = -2/rad_pts[j] *9/5

        mat_array[i+3,1] = j
        mat_array[i+3,2] = j+1 + n_rad_max
        mat_array[i+3,3] = -2 *3/5*df1_matsp[j,j+1]

    end

    for k = 1:tot_even_leg
        n_leg = 2*k
        for j = 3:n_rad_max-2

            i = 6*(j-3)+odd_entries_per_block+(k-1)*even_entries_per_row
            j0 = n_rad_max*(2*k-1)
            jm1 = j0 - n_rad_max
            jp1 = j0 + n_rad_max

            mat_array[i+1,1] = j +j0
            mat_array[i+1,2] = j-1+jm1
            mat_array[i+1,3] = 2 * (n_leg-1)/(2*n_leg-1) * df1_matsp[j,j-1]

            mat_array[i+2,1] = j +j0
            mat_array[i+2,2] = j+jm1
            mat_array[i+2,3] = -2/rad_pts[j] * (n_leg-1)^2/(2*n_leg-1) 

            mat_array[i+3,1] = j +j0
            mat_array[i+3,2] = j+1+jm1
            mat_array[i+3,3] = 2 * (n_leg-1)/(2*n_leg-1) * df1_matsp[j,j+1]

            mat_array[i+4,1] = j +j0
            mat_array[i+4,2] = j-1+jp1
            mat_array[i+4,3] = 2 * (n_leg+2)/(2*n_leg+3) * df1_matsp[j,j-1]

            mat_array[i+5,1] = j +j0
            mat_array[i+5,2] = j+jp1
            mat_array[i+5,3] = 2/rad_pts[j] * (n_leg+2)^2/(2*n_leg+3) 

            mat_array[i+6,1] = j +j0
            mat_array[i+6,2] = j+1+jp1
            mat_array[i+6,3] = 2 * (n_leg+2)/(2*n_leg+3) * df1_matsp[j,j+1]

        end

    end

    for k = 1:tot_even_leg-1
        n_leg = 2*k+1

        for j = 2:n_rad_max-1

            i = 6*(j-2)+odd_entries_per_block + tot_even_leg*even_entries_per_row + (k-1)*odd_entries_per_row
        
            j0 = n_rad_max*2*k
            jm1 = j0 - n_rad_max
            jp1 = j0 + n_rad_max

            mat_array[i+1,1] = j +j0
            mat_array[i+1,2] = j-1+jm1
            mat_array[i+1,3] = -2 * (n_leg-1)/(2*n_leg-1) * df1_matsp[j,j-1]

            mat_array[i+2,1] = j +j0
            mat_array[i+2,2] = j+jm1
            mat_array[i+2,3] = 2/rad_pts[j] * (n_leg-1)^2/(2*n_leg-1) 

            mat_array[i+3,1] = j +j0
            mat_array[i+3,2] = j+1+jm1
            mat_array[i+3,3] = -2 * (n_leg-1)/(2*n_leg-1) * df1_matsp[j,j+1]

            mat_array[i+4,1] = j +j0
            mat_array[i+4,2] = j-1+jp1
            mat_array[i+4,3] = -2 * (n_leg+2)/(2*n_leg+3) * df1_matsp[j,j-1]

            mat_array[i+5,1] = j +j0
            mat_array[i+5,2] = j+jp1
            mat_array[i+5,3] = -2/rad_pts[j] * (n_leg+2)^2/(2*n_leg+3) 

            mat_array[i+6,1] = j +j0
            mat_array[i+6,2] = j+1+jp1
            mat_array[i+6,3] = -2 * (n_leg+2)/(2*n_leg+3) * df1_matsp[j,j+1]
        end


    end


    n_leg = 2*tot_even_leg+1

    for j = 2:n_rad_max-1

        i = 3*(j-2)+odd_entries_per_block + tot_even_leg*even_entries_per_row + (tot_even_leg-1)*odd_entries_per_row
    
        j0 = n_rad_max*2*tot_even_leg
        jm1 = j0 - n_rad_max

        mat_array[i+1,1] = j +j0
        mat_array[i+1,2] = j-1+jm1
        mat_array[i+1,3] = -2 * (n_leg-1)/(2*n_leg-1) * df1_matsp[j,j-1]

        mat_array[i+2,1] = j +j0
        mat_array[i+2,2] = j+jm1
        mat_array[i+2,3] = 2/rad_pts[j] * (n_leg-1)^2/(2*n_leg-1) 

        mat_array[i+3,1] = j +j0
        mat_array[i+3,2] = j+1+jm1
        mat_array[i+3,3] = -2 * (n_leg-1)/(2*n_leg-1) * df1_matsp[j,j+1]

    end
    mat_array[even_entries_per_row*tot_even_leg+odd_entries_per_row*(tot_even_leg)+1,1] = (2*tot_even_leg+1)*n_rad_max
    mat_array[even_entries_per_row*tot_even_leg+odd_entries_per_row*(tot_even_leg)+1,2] = (2*tot_even_leg+1)*n_rad_max
    mat_array[even_entries_per_row*tot_even_leg+odd_entries_per_row*(tot_even_leg)+1,3] = 0.
    

    return mat_array
end


function calc_time_mat_array(n_rad_max::Int64,tot_even_leg::Int64,rad_pts::AbstractVector{Float64},df1_matsp::AbstractArray{Float64},df2_matsp::AbstractArray{Float64})

    tot_odd_entries = (tot_even_leg+1)*(n_rad_max-2)
    tot_even_entries = tot_even_leg*(n_rad_max-4)*3

    odd_entries_per_row = n_rad_max-2
    even_entries_per_row = 3*(n_rad_max-4)

    mat_array = Array{Float64,2}(undef,tot_even_entries+tot_odd_entries+1,3)

    for k = 1:tot_even_leg+1

        
        j0=2*(k-1)*n_rad_max
        for j = 2:n_rad_max-1

            i = j-1+(k-1)*odd_entries_per_row
            
            mat_array[i,1] = j + j0
            mat_array[i,2] = j + j0
            mat_array[i,3] = 1.
        end
    end

    for k = 1:tot_even_leg
        
        for j = 3:n_rad_max-2

            i = 3*(j-3)+tot_odd_entries + (k-1)*even_entries_per_row
            
            j0 = (2*k-1)*n_rad_max

            mat_array[i+1,1] = j+j0
            mat_array[i+1,2] = j-1+j0
            mat_array[i+1,3] = df2_matsp[j,j-1] + 2/rad_pts[j] * df1_matsp[j,j-1] 

            mat_array[i+2,1] = j+j0
            mat_array[i+2,2] = j+j0
            mat_array[i+2,3] = df2_matsp[j,j]  - 2*k*(2*k+1)/rad_pts[j]^2

            mat_array[i+3,1] = j+j0
            mat_array[i+3,2] = j+1+j0
            mat_array[i+3,3] = df2_matsp[j,j+1] + 2/rad_pts[j] * df1_matsp[j,j+1] 

        end

    end

    mat_array[tot_even_entries+tot_odd_entries+1,1] = (2*tot_even_leg+1)*n_rad_max
    mat_array[tot_even_entries+tot_odd_entries+1,2] = (2*tot_even_leg+1)*n_rad_max
    mat_array[tot_even_entries+tot_odd_entries+1,3] = 0.
    return mat_array

    

end

function calc_bc_mat_array(n_rad_max::Int64,tot_even_leg::Int64)

    tot_odd_entries = (tot_even_leg+1)*2
    tot_even_entries = tot_even_leg*4

    mat_array = Array{Float64,2}(undef,tot_even_entries+tot_odd_entries,3)

    mat_array[1,1] = 1
    mat_array[1,2] = 1
    mat_array[1,3] = 1.

    mat_array[2,1] = n_rad_max
    mat_array[2,2] = n_rad_max
    mat_array[2,3] = 1.

    for k = 1:tot_even_leg
        
        j1 = (2*k-1)*n_rad_max
        j0 = 2*k*n_rad_max
        i = 6*(k-1) + 2 

        mat_array[i+1,1] = j0 + 1
        mat_array[i+1,2] = j0 + 1
        mat_array[i+1,3] = 1.
        
        mat_array[i+2,1] = j0 + n_rad_max
        mat_array[i+2,2] = j0 + n_rad_max
        mat_array[i+2,3] = 1.

        mat_array[i+3,1] = j1 + 1
        mat_array[i+3,2] = j1 + 1
        mat_array[i+3,3] = 1.
        
        mat_array[i+4,1] = j1 + n_rad_max
        mat_array[i+4,2] = j1 + n_rad_max
        mat_array[i+4,3] = 1.

        mat_array[i+5,1] = j1 + 2
        mat_array[i+5,2] = j1 + 2
        mat_array[i+5,3] = 1.
        
        mat_array[i+6,1] = j1 + n_rad_max-1
        mat_array[i+6,2] = j1 + n_rad_max-1
        mat_array[i+6,3] = 1.




    end

    return mat_array


end

function mat_array_to_spmat(mat_array::AbstractArray{Float64})

    return sparse(mat_array[:,1],mat_array[:,2],mat_array[:,3])
end

function spatial_mats(n_rad_max::Int64,n_leg_max::Int64,rad_ratio::Float64)

    if iseven(n_leg_max)
        n_leg_max += 1 
    end

    r_i,r_o,rad_pts = calc_rad_pts(n_rad_max,rad_ratio)

    tot_even_leg::Int64 = (n_leg_max-1)÷2

    df1_matsp = df1_mat(n_rad_max,rad_pts)
    df2_matsp = df2_mat(n_rad_max,rad_pts)

    df3_matsp = *(df1_matsp,df2_matsp)
    df4_matsp = *(df2_matsp,df2_matsp)
    
    mat_array = calc_diff_mat_array(n_rad_max,tot_even_leg,rad_pts,df1_matsp,df2_matsp,df3_matsp,df4_matsp)

    diffusion_matsp = sparse(mat_array[:,1],mat_array[:,2],mat_array[:,3])
    
    mat_array = calc_cor_mat_array(n_rad_max,tot_even_leg,rad_pts,df1_matsp)

    coriolis_matsp = sparse(mat_array[:,1],mat_array[:,2],mat_array[:,3])

    mat_array = calc_time_mat_array(n_rad_max,tot_even_leg,rad_pts,df1_matsp,df2_matsp)

    time_matsp = sparse(mat_array[:,1],mat_array[:,2],mat_array[:,3])

    mat_array = calc_bc_mat_array(n_rad_max,tot_even_leg)

    bc_matsp = sparse(mat_array[:,1],mat_array[:,2],mat_array[:,3])


    return diffusion_matsp,coriolis_matsp,time_matsp,bc_matsp
end


function order_0_parametrized_system(n_rad_max::Int64,n_leg_max::Int64,rad_ratio::Float64,Ekman::Float64,freq::Float64)

    if iseven(n_leg_max)
        n_leg_max += 1 
    end

    diffusion_matsp,coriolis_matsp,time_matsp,bc_matsp = spatial_mats(n_rad_max,n_leg_max,rad_ratio)

    mat_G = 1im*freq*time_matsp+(1+0im)*(coriolis_matsp-Ekman*diffusion_matsp+bc_matsp)
    
    rhs = Vector{ComplexF64}(undef,n_leg_max*n_rad_max)
    
    rhs[:] .= 0
    rhs[n_rad_max] = -1
    

    vecs = mat_G \ rhs 
    
    writedlm("vec_jl_real.txt",real(vecs))
    writedlm("vec_jl_imag.txt",imag(vecs))

    return vecs,mat_G
end


function trapezoidal_int(vals::AbstractVector{Float64},n_x::Int64,x_pts::AbstractVector{Float64})

    int_result = 0.
    for k=1:n_x-1

        int_result += (vals[k]+vals[k+1])*(x_pts[k+1]-x_pts[k])

    end

    return int_result*0.5

end
