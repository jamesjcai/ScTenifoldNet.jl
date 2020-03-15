module scTenifoldNet

using Statistics, LinearAlgebra, Distributions, MultipleTesting, Random
import TSVD
import TensorToolbox

export tenrnet, manialn, drgenes

const NCOMP1,NCOMP2=3,5
const NLAYERS,NCELLS=10,500

# vecnorm(x) = x./norm.(x[:,i] for i in 1:size(x,2))'
vecnorm(x) = norm.(x[:,i] for i in 1:size(x,2))
function vecnorm!(x)
    for i in 1:size(x,2)
        x[:,i]=x[:,i]./norm(x[:,i])
    end    
end


function pcnet(X::AbstractArray{T,N}, p::Int=3;
              scalein=false, scaleout=false, symmout=false) where N where T<:Real
    if scalein
        σ=std(X,dims=1)
        σ(σ.==0).=1
        X=(X.-mean(X,dims=1))./σ
    end
    ℊ=size(X,2)
    A=1.0 .-Matrix(I,ℊ,ℊ)
    Threads.@threads for k in 1:ℊ
        y=X[:,k]
        𝒳=X[:,1:end.≠k]
        ϕ=TSVD.tsvd(𝒳,p)[3];
        s=𝒳*ϕ
        s ./= (vecnorm(s).^2)'
        b=sum(y.*s,dims=1)
        𝒷=ϕ*b'
        @inbounds A[k,A[k,:].==1.0]=𝒷
    end
    if symmout
        A=0.5*(A+A')
    end
    if scaleout
        A=A./maximum(abs.(A))
    end
    return convert(Array{Float16,2},A)
  end

function tensordecomp(Λ, p::Int64=5)
    𝒯=TensorToolbox.cp_als(Λ,p)
    𝕏=TensorToolbox.full(𝒯)
    A=mean(𝕏[:,:,i] for i=1:size(𝕏,3))
    # A ./=maximum(abs.(A))
    # A=round.(A; digits=5)
    return A
end

function manialn(X,Y)
    μ,dim=0.9,30
    n1,n2=size(X,1),size(Y,1);
    W₁,W₂=X.+1,Y.+1
    ℐ=Matrix(I,n1,n2);
    μ = μ*(sum(W₁)+sum(W₂)/(2*sum(ℐ)));
    𝕎 = [W₁ μ*ℐ; μ*ℐ' W₂];
    L=diagm(vec(sum(abs.(𝕎),dims=1))).-𝕎;
    λ,V = eigen(L);
    i=real(λ).>=1e-8;
    V=real(V[:,i]);
    V=V[:,1:dim];    
    aln1=V[1:n1,:];
    aln2=V[n1+1:end,:];
    d = norm.((aln1.-aln2)[i,:] for i = 1:n1)
    # _,idx=findmax(dd)
    return d, aln1, aln2
end

function drgenes(d)
    d²=d.^2
    FC=d²./mean(d²)
    pVals = ccdf.(Chisq(1),FC)
    pAdjusted = MultipleTesting.adjust(pVals, BenjaminiHochberg())
    return FC,pVals,pAdjusted
end

function tenrnet(X::AbstractMatrix{T}; donorm::Bool=true) where T
    ℊ,𝒸=size(X)
    if donorm
        lbsz=sum(X,dims=1)
        X=(X./lbsz)*median(lbsz)
        # X=(X./lbsz)*1e4
    end    
    A=zeros(Float16, ℊ, ℊ, NLAYERS)
    for k=1:NLAYERS
        println("network ... $k")
        𝕩=X[:,randperm(𝒸)][:,1:NCELLS]
        @time @inbounds A[:,:,k]=pcnet(𝕩',NCOMP1)
    end
    Z=tensordecomp(A,NCOMP2)
    return Z
end

end # module
