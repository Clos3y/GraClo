
GraClo := module()

export IndependentComponents,NumberOfIndependentComponents,Equations,SYM,ASYM;
option package;
SYM := proc(condition::list,tenss);

local perms := [],p,l,k,base,ind;

base := cat(seq(parse(convert(tenss,string)[i]),i=1...StringTools:-Search("[",convert(tenss,string)) - 1)):
ind := op(tenss):

for p in combinat:-permute(condition) do:

p := convert(p,Array):

for l,k in ind do:

   if member(k,condition) then else 

   p := ArrayTools:-Insert(p,l,ind[l]);
   
   end if;

od:

perms := [op(perms),convert(p,list)];

od:

return (1/nops(perms))*add(base[op(perms[i])],i=1...nops(perms));

end proc:

ASYM := proc(condition::list,tenss);

local perms := [],p,l,k,base,ind;

base := cat(seq(parse(convert(tenss,string)[i]),i=1...StringTools:-Search("[",convert(tenss,string)) - 1)):
ind := op(tenss):

for p in combinat:-permute(condition) do:

p := convert(p,Array):

for l,k in ind do:

   if member(k,condition) then else 

   p := ArrayTools:-Insert(p,l,ind[l]);
   
   end if;

od:

perms := [op(perms),convert(p,list)];

od:

return (1/nops(perms))*add(`if`(i mod 4 < 2,base[op(perms[i])],-1*base[op(perms[i])]),i=1...nops(perms));

end proc:

IndependentComponents := proc(tenss,basisEquations::list,dim::posint)
local sol:=[],dummyList:=[],base,ind,rank,k,perm,i,solutions,depList,subCond,finalEquationList;

base := cat(seq(parse(convert(tenss,string)[i]),i=1...StringTools:-Search("[",convert(tenss,string)) - 1)):
ind := op(tenss):
rank := nops([ind]):

for k from 0 to (dim^rank - 1) do:


 perm[k] := op(ListTools:-Reverse(convert(k,'base',dim))):


 if numelems([perm[k]]) < rank then


  for i from 1 to (rank - numelems([perm[k]])) do:


  perm[k] := 0,perm[k]


  od:


 else end if:


 dummyList := [op(dummyList),[perm[k]]]:


od: 

solutions := []:
depList := []:
if nops(basisEquations) = 0 then solutions := [seq(base[op(op(dummyList)[m])],m=1...numelems(dummyList))] else
subCond := [seq(ind =~ dummyList[i],i=1..numelems(dummyList))]:
finalEquationList := ListTools:-MakeUnique(map(op,[seq(`if`(subs(subCond[i],basisEquations) = [0=0],NULL,subs(subCond[i],basisEquations)),i=1...numelems(subCond))])):
for i in solve(finalEquationList,maxsols=infinity) do:
if lhs(i) = rhs(i) then solutions:=[op(solutions),lhs(i)] else depList :=[op(depList),lhs(i)] end if;
od:
end if:

return solutions;
end proc;

NumberOfIndependentComponents := proc(tenss,basisEquations::list,dim::posint)
local sol:=[],dummyList:=[],base,ind,rank,k,perm,i,solutions,depList,subCond,finalEquationList;

base := cat(seq(parse(convert(tenss,string)[i]),i=1...StringTools:-Search("[",convert(tenss,string)) - 1)):
ind := op(tenss):
rank := nops([ind]):

for k from 0 to (dim^rank - 1) do:


 perm[k] := op(ListTools:-Reverse(convert(k,'base',dim))):


 if numelems([perm[k]]) < rank then


  for i from 1 to (rank - numelems([perm[k]])) do:


  perm[k] := 0,perm[k]


  od:


 else end if:


 dummyList := [op(dummyList),[perm[k]]]:


od: 

solutions := []:
depList := []:
if nops(basisEquations) = 0 then solutions := [seq(base[op(op(dummyList)[m])],m=1...numelems(dummyList))] else
subCond := [seq(ind =~ dummyList[i],i=1..numelems(dummyList))]:
finalEquationList := ListTools:-MakeUnique(map(op,[seq(`if`(subs(subCond[i],basisEquations) = [0=0],NULL,subs(subCond[i],basisEquations)),i=1...numelems(subCond))])):
for i in solve(finalEquationList,maxsols=infinity) do:
if lhs(i) = rhs(i) then solutions:=[op(solutions),lhs(i)] else depList :=[op(depList),lhs(i)] end if;
od:
end if:

return nops(solutions);
end proc;

Equations := proc(tenss,basisEquations::list,dim::posint)
local sol:=[],dummyList:=[],base,ind,rank,k,perm,i,solutions,depList,subCond,finalEquationList,deps;

base := cat(seq(parse(convert(tenss,string)[i]),i=1...StringTools:-Search("[",convert(tenss,string)) - 1)):
ind := op(tenss):
rank := nops([ind]):

for k from 0 to (dim^rank - 1) do:


 perm[k] := op(ListTools:-Reverse(convert(k,'base',dim))):


 if numelems([perm[k]]) < rank then


  for i from 1 to (rank - numelems([perm[k]])) do:


  perm[k] := 0,perm[k]


  od:


 else end if:


 dummyList := [op(dummyList),[perm[k]]]:


od: 

solutions := []:
depList := []:
if nops(basisEquations) = 0 then solutions := [seq(base[op(op(dummyList)[m])],m=1...numelems(dummyList))] else
subCond := [seq(ind =~ dummyList[i],i=1..numelems(dummyList))]:
finalEquationList := ListTools:-MakeUnique(map(op,[seq(`if`(subs(subCond[i],basisEquations) = [0=0],NULL,subs(subCond[i],basisEquations)),i=1...numelems(subCond))])):
for i in solve(finalEquationList,maxsols=infinity) do:
if lhs(i) = rhs(i) then solutions:=[op(solutions),lhs(i)] else depList :=[op(depList),lhs(i)] end if;
od:
end if:

deps := op(solve(finalEquationList,depList)):
depList := [seq(`if`(rhs(deps[i])=0,NULL,deps[i]),i=1...nops(deps))]:
return depList;

end proc;

end module:
