using DelimitedFiles

D21=Dict{Int64,Char}(1 =>'-',2 =>'A',3 =>'C',4 =>'D',5 =>'E',6 =>'F',7 =>'G',8 =>'H',9 =>'I',10 =>'K',11 =>'L',12 =>'M',13 =>'N',14=>'P',15=>'Q',16 =>'R',17=>'S',18 =>'T',19 =>'V',20 =>'W',21 =>'Y')
Di21=Dict{Char,Int64}('-'=>1,'X'=>1,'B'=>1,'J'=>1,'U'=>1,'Z'=>1,'O'=>1,'A'=>2,'C'=>3,'D'=>4,'E'=>5,'F'=>6,'G'=>7,'H'=>8,'I'=>9,'K'=>10,'L'=>11,'M'=>12,'N'=>13,'P'=>14,'Q'=>15,'R'=>16,'S'=>17,'T'=>18,'V'=>19,'W'=>20,'Y'=>21)

function readmsanum(infile::String ; format=1, header=false)
	Y = Array{Float64,2}(undef,0,0)
	try 
		if header
			Y = readdlm(infile, Int64, skipstart=1)
		else
			Y = readdlm(infile, Int64)
		end	
	catch err
		println("inputoutput.jl - readmsanum: readdlm failed on file $infile. The alignment may not be of the correct format.")
		error(err)
	end

	if format==0
		Y .+= 1
	elseif format==1
		if findmin(Y)[1] == 0
			error("inputoutput.jl - readmsanum: file $infile contains a 0 value. `format` should be set to 0.")
		end
	end

	return Y
end



function transform_seq(seq::Any,D::Dict{Char,Int64})
    q=length(seq)
    let_seq=Array{Int64}(undef,q)
    for i in 1:q
       let_seq[i]=D[seq[i]]
    end    
  return let_seq
end 


#read the MSA as fasta file and transform to a MSA as a matrix of numbers in {1,...,21}.
function translate_fasta_to_num_matrix(msa_fasta::String)
    pfam_sa = readlines(msa_fasta);
    M=Int(length(pfam_sa)/2)
    L=length(pfam_sa[2])
	let_msa=Array{Int64,2}(undef,M,L);
    index=1
	for i in 1:length(pfam_sa)
        if i%2==0
         let_msa[index,:].=transform_seq(pfam_sa[i],Di21)
            index+=1
		end
	    
	end
return let_msa
end


function export_fasta_file(outfile::String,msa::Array{Char,2})
(N,L)=size(msa)
fname=outfile
open(fname,"w") do file
      for i in 1:N
       write(file,">sequence_$i")
       write(file,"\n")
       write(file,String(msa[i,:]))
       write(file,"\n") 
       end

end

end

#transform MSA from num to letters
function transform_MSA_fasta(msa::Array{Int64,2})
	(N,L)=size(msa)
	let_msa=Array{Char}(undef,N,L);
	for i in 1:N
		for j in 1:L
		 let_msa[i,j]=D21[msa[i,j]]   
        end
	    
	end
return let_msa
end
