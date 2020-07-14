module Null_models_I_and_II
using LinearAlgebra


export sample_from_Null_model_I, sample_from_Null_model_II, PairwiseHammingDist,translate_fasta_to_num_matrix, transform_MSA_fasta, export_fasta_file
include("utils.jl")

##Dirac delta 
function deltak(i::Int64,j::Int64)
  i==j ? 1.0 : 0.0
    end

###compute hamming distance matrix D
 function PairwiseHammingDist(sample::Array{Int64,2})
    Y=sample[1:end,1:end]
    Y = Y'
    (N,M) = size(Y)
    hav = zeros(M,M)
    
  @inbounds  for m in 1:M
               for n in (m+1):M
                hd = 0.
                for i = 1:N
                 hd += Y[i,m]!=Y[i,n]
                 end
               hd = hd / N
            
              hav[m,n] = hav[n,m]=hd
               end
              end
    return hav
end

function PairwiseHammingDist(sample::Any)
   if typeof(sample)==Array{Char,2}
     (M,N)=size(sample)
	 Y=Array{Int}(undef,M,N);
	 for i in 1:M
		for j in 1:N
         Y[i,j]=Di21[sample[i,j]]   
        end
     end
    elseif typeof(sample)==Array{Int64,2}
          Y=sample[1:end,1:end]
    elseif typeof(sample)==String
        sample=translate_fasta_to_num_matrix(sample)
        Y=sample[1:end,1:end]
    else        
     error("Null_models_I_and_II.jl - PairwiseHammingDist: input parameter format is wrong.\n") 
   end              
    Y = Y'
    (N,M) = size(Y)
    hav = zeros(M,M)
    
  @inbounds  for m in 1:M
               for n in (m+1):M
                hd = 0.
                for i = 1:N
                 hd += Y[i,m]!=Y[i,n]
                 end
               hd = hd / N
            
              hav[m,n] = hav[n,m]=hd
               end
              end
    return hav
end

#given column k and two rows m,n an exchange of the align entries a_km and a_kn is attempted.
function swap_align!(align::Array{Int64,2},m1::Int64,m2::Int64,k::Int64)
    alpha=align[m1,k]
    bet=align[m2,k]
    align[m1,k]=bet
    align[m2,k]=alpha
end

### cost function
function Energy(dtarget::Array{Float64,2},dnew::Array{Float64,2})
 return sum((dtarget .- dnew) .^2)
end 

## frobenious distance between matrices  
function compute_l2_error(d1::Array{Float64,2},d2::Array{Float64,2})
     
     return sqrt(sum((d1 .- d2).^2))/size(d1)[1]
end 

# return delta matrix which store the increment in the hamming distance entries
function delta_d(align::Array{Int64,2},delta::Array{Float64,2},m1::Int64,m2::Int64,k::Int64)
    (M,L)=size(align)
    #delta .=0
 @inbounds @simd for j in 1:M
        if j!=m1&&j!=m2
        delta[j,m1]=-(deltak(align[m2,k],align[j,k])-deltak(align[m1,k],align[j,k]))/L
        delta[j,m2]=-delta[j,m1]
        delta[m1,j]=delta[j,m1]
        delta[m2,j]=delta[j,m2]    
            end 
      
        
        end 
    
    
   return delta 
end

function delta_d!(align::Array{Int64,2},delta::Array{Float64,2},m1::Int64,m2::Int64,k::Int64)
    (M,L)=size(align)
 @fastmath @inbounds @simd for j in 1:M
        if j!=m1&&j!=m2
        delta[j,m1]=-1.0*(deltak(align[m2,k],align[j,k])-deltak(align[m1,k],align[j,k]))/L
        delta[j,m2]=-1.0*delta[j,m1]
        delta[m1,j]=delta[j,m1]
        delta[m2,j]=delta[j,m2]    
            end 
     
        end 
    
end
  
## change produced at d_current as consequence of swap entries m1 and m2
 function change_d!(d_current::Array{Float64,2},delta::Array{Float64,2},m1::Int64,m2::Int64)
        M=size(d_current)[1]

    @fastmath @inbounds @simd for j in 1:M
                     if j!=m1&&j!=m2
                     d_current[j,m1]=d_current[j,m1]+delta[j,m1]
                     d_current[j,m2]=d_current[j,m2]+delta[j,m2]
                     d_current[m1,j]=d_current[j,m1]
                     d_current[m2,j]=d_current[j,m2]
                     end
                    end

        end

##Frobenious norm betwenn hamming distances matrices at times t+1 (dinic) respect to dtarget
## given an alignment permutation
 function dE(dtarget::Array{Float64,2},dinic::Array{Float64,2},delta::Array{Float64,2},m1::Int64,m2::Int64)
  output=0.0    
  s=size(dtarget)[1]
  @fastmath @inbounds for j in 1:s 
  if j!=m1&&j!=m2
  output+=2*(dinic[j,m1]-dtarget[j,m1])*delta[j,m1]+delta[j,m1]^2+2*(dinic[j,m2]-dtarget[j,m2])*delta[j,m2]+delta[j,m2]^2                 
              end
       end
    return 2.0*output
 end

   function do_swap(align::Array{Int64,2},delta::Array{Float64,2},d_target::Array{Float64,2},d_current::Array{Float64,2},temp::Float64,energ::Float64)
    
    (M,L)=size(align)
    m1=rand(1:M)
    m2=rand(1:M)
    k=rand(1:L)
    de=0.0
    while align[m1,k]==align[m2,k]
           m1=rand(1:M)
           m2=rand(1:M)
           k=rand(1:L)  
    end 


     delta_d!(align,delta,m1,m2,k)###
     de=dE(d_target,d_current,delta,m1,m2)####
     if de<0.0 || rand()< exp(-de/temp)
            swap_align!(align,m1,m2,k)
            energ=energ+de
            change_d!(d_current,delta,m1,m2)###
               return 1,energ,d_current
           else
            return 0,energ,d_current
           end

 
    end      

function update_param_SA(alignment::Array{Int64,2},alignment_best::Array{Int64,2},n_it_tot::Int64,n_it_current::Int64,flag::Int64,acceptance_rate::Float64, energy_min ::Float64,energy_min_prev ::Float64,energy::Float64,T_init::Float64,T_min::Float64,n_it_max::Int64,T_factor_slow::Float64,T_factor_fast::Float64)
 (M,L)=size(alignment)
 stop_condition = 0
 n_it_max =n_it_max#N_IT_MAX;
 flag_max = 100;   
 n_swap_min = convert(Int64,round(M*L/5))
 acceptance_threshold_low = 0.0001 
 acceptance_threshold_high = 0.95
 T=T_init
 T_min = T_min#2/(10*M)
 #T_factor_slow =T_factor_slow#T_FACTOR_SLOW, 
 #T_factor_fast =T_factor_slow#T_FACTOR_FAST;

 #Update number of iterations for next MCMC
    if acceptance_rate!=0
      n_it_current = convert(Int64,round(n_swap_min / acceptance_rate))
    else
        n_it_current = 2*n_it_current
    end    
        
    if n_it_current > M*L
        n_it_current = M*L
    end
##If a better solution has been found and acceptance rate is high enough, decrease temperature
    if energy_min_prev > energy_min && acceptance_rate > acceptance_threshold_low
        
#linear cooling scheme with two slopes
        if acceptance_rate > acceptance_threshold_high
        T = T*T_factor_fast;
        else
        T = T*T_factor_slow
        end    
 #      
        flag -=1;
        if flag<0
        flag=0
        end    
    ## Else, give it another chance at this temperature, starting from the best solution
    else
        alignment .= alignment_best
        energy = energy_min;
        flag+=1
    end
    ## Test various stop conditions
    if T<T_min 
    stop_condition=1
     println("SA terminated -- Temperature reached its minimum")
    elseif n_it_tot>n_it_max
        stop_condition = 1 
     println("SA terminated -- Maximum number of iterations")
        
    elseif flag>flag_max
        stop_condition=1 
    println("SA terminated - Could not find better solution OR acceptance rate too low $(flag).")
        end
    return (stop_condition,alignment,n_it_current,flag,energy,T);
end

function shuffle_alignment(alignment::Array{Int64,2},T::Float64,d_curr::Array{Float64,2}, d_target::Array{Float64,2})

 (M,L)=size(alignment)
 n_it_tot=0
 n_it_curr = 10000 
 n_it_max=2*M*L
 stop_condition = 0
 stop_condition_max = 100
 err_l2 = compute_l2_error(d_curr, d_target)
 d_pos=zeros(M,M)

    while stop_condition<stop_condition_max
        n_swap_done = 0
        for i in 1:n_it_curr
         (add,energy,d_curr) = do_swap(alignment,d_pos,d_target,d_curr,T,err_l2)   
             n_swap_done = n_swap_done + add;
        end
        acceptance_rate = n_swap_done/(n_it_curr-1)
        n_it_tot = n_it_tot + n_it_curr;
        err_l2 = compute_l2_error(d_curr,d_target);
        if acceptance_rate < 0.99
            T = 2*T
        else
            T = T / 1.1
        end
        println("n_it_tot:$(n_it_tot)")
        if n_it_tot>n_it_max
            stop_condition=stop_condition_max
        end  
    end 
    println("$(n_it_tot) iterations for shuffling\n")
    return T,d_curr
end 

### run simmulating annealing with a linear cooling scheme with two slopes
function run_SA(align::Array{Int64,2},T_init::Float64,T_min::Float64,T_factor_slow::Float64,T_factor_fast::Float64,shuffle::Int64,num_iter_max::Int64,d_curr::Array{Float64,2},d_target::Array{Float64,2})
    
    (M,L)=size(align)
    T = T_init
    energy=0.0 
    energy_min=0.0
    energy_min_prev=100.0
    acceptance_rate = 0.0; 
    n_it_max = num_iter_max
    n_it_current = 1
    n_it_tot = 0
    stop_condition = 0 
    flag=0 
    ## Initial state
    println("runMC\n")
    alignment_best=align[1:end,1:end]
    err_l2 = compute_l2_error(d_curr,d_target)
    err_l2_min = err_l2
    energy = Energy(d_curr, d_target)

    energy_min = energy
    #println("Initial energy : $(energy) -- Initial l2 error: $(err_l2)")

    d_i=d_curr[1:end,1:end]


    if shuffle!=0
    
        println("Shuffling system - T = $(T) ...")
        (T,d_curr)=shuffle_alignment(align,T,d_i,d_target)
        err_l2 = compute_l2_error(d_curr, d_target)
        err_l2_min = err_l2
        energy = Energy(d_curr, d_target)
        energy_min = energy;
        println("Done.");
        println("After shuffling : energy = $(energy) -- l2error=$(err_l2) -- T =$(T)");
    end

    ## Sampling
    println("\nSampling ... \n")
        d_pos=zeros(M,M)

    while stop_condition==0  
     #Calculation 
    println("n_it_current: $(n_it_current)")
    println(" T : $(T)");
    println("energy : $(energy)");
    #println("energy 2 : $(Energy(d_target,d_curr))")
    println("acceptance_rate : $(acceptance_rate)");
      n_swap_done = 0;
           for i in 1:n_it_current
           	
         (add,energy,d_curr) = do_swap(align,d_pos,d_target,d_curr,T,energy)  

        n_swap_done = n_swap_done + add;
             if energy < energy_min
             energy_min_prev = energy_min 
             energy_min = energy
             acceptance_rate = n_swap_done/i;
                if acceptance_rate < 0.001 
             alignment_best .= align
                end    
             end
          end
        
        acceptance_rate =  n_swap_done/n_it_current
        n_it_tot = n_it_tot + n_it_current
        err_l2 = compute_l2_error(d_curr,d_target)
        
        if err_l2<err_l2_min
          err_l2_min = err_l2
        end
        ## Just to have alignment_best initialized
        if acceptance_rate > 0.001 
            alignment_best .= align
         end 
  println("n_it_tot:$(n_it_tot)")

 ## Updating parameters
 (stop_condition,align,n_it_current,flag,energy,T)=update_param_SA(align,alignment_best,n_it_tot,n_it_current,flag,acceptance_rate, energy_min,energy_min_prev,energy,T,T_min,n_it_max,T_factor_slow,T_factor_fast)   
 
    end

    println(" \ndone\n");
    println("n_it_tot:$(n_it_tot)")
    println("Final l2 error = $(err_l2_min)")

  return align

end

####### given a MSA in numeric format this function provide a randomized MSA according to Null-model I
function sample_from_Null_model_I(infile::String;outfile="",shuffle_temp=10.0)
  msa=translate_fasta_to_num_matrix(infile)	
  println("Shuffling alignment")
  
  dtarget=PairwiseHammingDist(msa);
  shuffle_alignment(msa,shuffle_temp,dtarget,dtarget)
  msa_fasta=transform_MSA_fasta(msa)
  if !isempty(outfile)
  export_fasta_file(outfile,msa_fasta)
  println("randomized MSA has been exported with fasta format!!")
  end

 return msa_fasta
end



####### given a MSA in numeric format this function provide a randomized MSA according to Null-model II
function sample_from_Null_model_II(infile::String;outfile="",shuffle=1,shuffle_temp=10.0,T_factor_slow=0.8,T_factor_fast=0.1,min_temp="default",num_iter_max=20000000)
  msa=translate_fasta_to_num_matrix(infile)
  dtarget=PairwiseHammingDist(msa);
  (M,L)=size(msa)
  if min_temp=="default"
  	t_min=2/(10*M)
  else
  	t_min=min_temp
  end
  println("Parameters:")
  if shuffle==1
     println("Simulating annealing start from shuffled alignment")
  end 
  println("minimal temperature  T_min=$(t_min)")
  println("slow cooling rate = $(T_factor_slow), fast cooling rate = $(T_factor_fast)")
  println("maximum number of iterations $(num_iter_max)")



  result=run_SA(msa,shuffle_temp,t_min,T_factor_slow,T_factor_fast,shuffle,num_iter_max,dtarget,dtarget)
  msa_fasta=transform_MSA_fasta(result)
  if !isempty(outfile)
  export_fasta_file(outfile,msa_fasta)
  println("randomized MSA has been exported with fasta format!!")
  end
 return msa_fasta
end



end
