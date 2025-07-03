
restart; usetimes:=false:
with(ListTools): with(LinearAlgebra): with(Threads): with(Threads[Task]): with(ArrayTools): with(plots): with(ListTools):
# The goal of this maple file is to solve the generalised eigenvalue problem for arbitrary multiple backflow operators, as detailed in "Repeated quantum backflow and overflow" (https://arxiv.org/abs/2505.13184) using tiral vectorspsi__n(q) = q^(n + delta)*e^(-Aq);where abs(delta) < 1/2, 0 < A;are parameters to be fixed. 
# 
# The parameters below mean as follows:
# N; - Largest trial vector used in the sense that the largest backflow and overflow spectral points are found in the subspace of L^2;functions spanned by `(&psi;__n)__0&le;n&le;N`;.
# usetimes - A boolean variable that is used to decide whether to use custom times or the default times outlined below. 
# M; - Number of backflow periods under investigation
# delta; - Shift of the exponent of the trial vectors as defined above
# A; - Exponential multiplier of the trial vectors as defined above
# 
# ***By default the Maple code works with the assumption that the backflow times under investigation are given by t__k = (2*k - 2*M + 1)/2;for k = 0 .. 2*M - 1.;Particular times can be specified using the times; variable.***
# dps;- Digits of the file, by default this is set to N;.
workingdir:=currentdir();
currentdir(workingdir):
if not(N::integer) then N:=100; fi;
if not(assigned('delta')) then delta:=0; fi;
if not(assigned('A')) then A:=1; fi;
if not(Nstart::integer) then Nstart:=1; fi;
if not(Nstop::integer) then Nstop:=100; fi;
if not(assigned(dps)) then dps:=N; fi;
if A<1 then
tA:=cat("0",A):
else
tA:=A:
fi;
if delta=0 then
tdelta:="0.0":
elif delta<1 and delta>0 then
tdelta:=cat("0",delta):
elif delta<0 then
tdelta:=cat("-0",abs(delta));
else
tdelta:=delta:
fi;
if not(M::integer) then M:=3; 
if not(usetimes) then MatrixElement_string:=cat("auxvals/MB_N0_",N,"MAX",N,"_M",M,"_delta",tdelta,"_a",tA,"_DPS",dps,".mpl"); fi;
fi;
if not(assigned(times)) and usetimes then times:="-0.901,-0.346,-0.099,0.099,0.346,0.901"; 
MatrixElement_string:=cat("auxvals/MB_N0_",N,"MAX",N,"_t",times,"_delta",delta,"_a",tA,"_DPS",dps,".mpl");fi;
dps:=N: Digits:=dps+40:
data_python:='data_python':
read MatrixElement_string:
Digits:=N+20:
shiftpar:=0.039:
split:=proc(n,k) local r,s;
s:=[0,seq(isqrt(floor(n*(n+1)*r/k)),r=1..k-1),n];
return([seq(s[r]+1..s[r+1],r=1..k)]);
end proc;
gram_element:=(n,m)->data_python['gram_element',n,m];
matrix_element:=(n,m)->piecewise(n=m,Re(data_python['matrix_element',n,m]),data_python['matrix_element',n,m]);
B:=Matrix(Nstop,Nstop,shape=hermitian,storage=triangular[lower]):
P:=Matrix(Nstop,Nstop,shape=symmetric,storage=triangular[lower]):
for ii from 1 to Nstop do:
for jj from 1 to ii do:
B[ii,jj]:=matrix_element(ii-1,jj-1):
P[ii,jj]:=gram_element(ii-1,jj-1):
od: od:
printf("Matrices loaded in.\n");
find_kth:=proc(L,k);
LL:=convert(L,list):
if k>nops(LL) or k<-nops(LL) or k=0 or not(k::integer) then
printf("Parameter k must be an integer with absolute value between 1 and the size of the list\n");
return 0
fi: 
LL_sorted:=sort(LL,numeric):
index:=Search(LL_sorted[k],L):
return LL[index],index;
end proc;
gen_eigs:=proc(B,P) local Btrunc,Ptrunc,maxests,minests; option remember;
maxests:=[]: minests:=[]: max2ests:=[]: min2ests:=[]: 
for ii from Nstart to Nstop do:
Btrunc:=B[1..ii,1..ii]: Ptrunc:=P[1..ii,1..ii]:
eigs:=Eigenvalues(Btrunc,Ptrunc):
maxeig:=max(Re(eigs)): max2eig:=find_kth(Re(eigs),-2)[1]:
mineig:=min(Re(eigs)): min2eig:=find_kth(Re(eigs),2)[1]:
printf("maxest %d: %f\n", ii, evalf[5](maxeig));
printf("max2est %d %f\n", ii, evalf[8](max2eig));
printf("minest %d: %f\n", ii, evalf[5](mineig));
printf("min2est %d %f\n", ii, evalf[8](min2eig));
maxests:=[op(maxests), maxeig]: minests:=[op(minests), mineig]:
max2ests:=[op(max2ests), max2eig]: min2ests:=[op(min2ests), min2eig]:
od:
return maxests, max2ests, minests, min2ests;
end proc;
maxests, max2ests, minests, min2ests:=gen_eigs(B, P):
if usetimes then
save maxests, cat("maxests/maxestsN",Nstart,"_",Nstop,"MAX",N,"M",M,"t",times,"d",tdelta,"A",A,".m");
save minests, cat("maxests/minestsN",Nstart,"_",Nstop,"MAX",N,"M",M,"t",times,"d",tdelta,"A",A,".m");
save max2ests, cat("maxests/max2estsN",Nstart,"_",Nstop,"MAX",N,"M","t",times,"d",tdelta,"A",A,".m");
save min2ests, cat("maxests/min2estsN",Nstart,"_",Nstop,"MAX",N,"M","t",times,"d",tdelta,"A",A,".m");
else 
save maxests, cat("maxests/maxestsN",Nstart,"_",Nstop,"MAX",N,"M",M,"d",tdelta,"A",A,".m");
save minests, cat("maxests/minestsN",Nstart,"_",Nstop,"MAX",N,"M",M,"d",tdelta,"A",A,".m");
save max2ests, cat("maxests/max2estsN",Nstart,"_",Nstop,"MAX",N,"M","d",tdelta,"A",A,".m");
save min2ests, cat("maxests/min2estsN",Nstart,"_",Nstop,"MAX",N,"M","d",tdelta,"A",A,".m");
listplot(maxests);
maxests[-1];
listplot(minests);
minests[-1];
listplot(min2ests);
min2ests[-1];
listplot(max2ests);
