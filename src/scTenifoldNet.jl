module scTenifoldNet

using Statistics, LinearAlgebra, Arpack, TensorToolbox, Random
export pcnet

function pcnet(X)
    n=size(X,2)
    A=1.0 .-Matrix(I,n,n)
    for k in 1:n        
        y=X[:,k]
        𝒳=X[:,1:end.≠k]
        _,ϕ=Arpack.eigs(𝒳'𝒳,nev=3,which=:LM)
        s=𝒳*ϕ
        s ./=((norm.(s[:,i] for i=1:size(s,2))).^2)'
        b=sum(y.*s, dims=1)
        𝒷=ϕ*b'
        A[k,A[k,:].==1.0]=𝒷
    end
    return A
end

function tensordecomp(X)
    a=TensorToolbox.cp_als(X,5)
    Z=full(a)
    Z1=mean(Z[:,:,i] for i=1:size(Z,3))
    Z1=Z1./maximum(abs.(Z1));
    # Z1=round.(Z1; digits=5)
    return Z1
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
    dd = norm.((aln1.-aln2)[i,:] for i = 1:n1)
    # _,idx=findmax(dd)
    return sortperm(-dd)    
end

function dddtenifold(X)
    lbsz=sum(X,dims=1)
    X=(X./lbsz)*median(lbsz)
    ℊ,𝒸=size(X)
    A=zeros(Float64, ℊ, ℊ, 10)
    for k=1:10
        println("network ... $k")
        𝕩=X[:,randperm(𝒸)][:,1:500]
        A[:,:,k]=pcnetworkhelp(𝕩')
    end
    Z=tensordecomp(A)
    return Z
end

end # module
