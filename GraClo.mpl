GraClo := module()

   description "A module for finding the independent components and corresponding dependent components (with dependence) given a set of equations.";

   #v2.1

   # These are the functions accessible outside of the module
   export IndependentComponents,NumberOfIndependentComponents,Equations,SYM,ASYM;

   #This packages the module
   option package;

   permuter := proc(condition::list,tensor)::Array:

      options remember,threadsafe;

      uses combinat,ArrayTools;

      #Initialises the local variables
      local ind :: exprseq := op(tensor);

      local permutations :: Array := Array(map(convert,permute(condition),Array),datatype=Array);

      #If the index appears in the tensor, then it is ignored. Else, it is added in the correct, previous position
      for local l,k in ind do:
       map[inplace](i -> `if`(member(k,condition),i,Insert(i,l,k)),permutations);
      od:

      return permutations;

   end proc;

   SYM := proc(condition::list,tensor)::`+`;

      description "Computes the symmetric sum over the specified indices";

      options remember,threadsafe;

      uses combinat;

      local base :: symbol := map2(op,0,tensor);
     
      local permutations := permuter(condition,tensor);

      local n_perm :: posint := numbperm(condition);
         #Returns the symmetric sum
       return (1/n_perm)*add(base[op(convert(permutations[i],list))],i=1...n_perm);

   end proc:

   ASYM := proc(condition::list,tensor)::`+`;

      description "Computes the antisymmetric sum over the specified indices";

      options remember,threadsafe;

      uses combinat,GroupTheory;

      local base :: symbol := map2(op,0,tensor);
     
      local permutations :: Array := permuter(condition,tensor);

      local n_perm :: posint := numbperm(condition);

      local signs :: list := map(PermParity,map(Perm,permute(nops(tensor))));

     return (1/n_perm)*add(signs[i]*base[op(convert(permutations[i],list))],i=1...n_perm);

   end proc:

   IndependentComponents := proc(tenss,basisEquations::list,dim::posint,startdim::integer:=0)::list;

      description "Find the independent components of a tensor given some equations defining its behaviour, and a dimension. One may optionally change the starting index";

      option remember,threadsafe;

      uses ListTools;

      #Initialises the local variables
      local sol:=[],dummyList:=[],base,ind,rank,k,perm,i,solutions,depList,subCond,finalEquationList;

      #These extract the base, indices and rank from the tensor. I can't yet find a better way to do this
      base := map2(op,0,tenss):
      ind := op(tenss):
      rank := nops([ind]):

      #Creates the permuted list of all numbered elements (i.e., [0,0,0],[0,0,1],[0,0,2]...)
      dummyList := [seq([seq(iquo(i,dim^r) mod dim^r,r=(rank-1)..1,-1),irem(i,dim)],i=0..(dim^(rank) - 1))];

      #Solves the equations for the independent components
      solutions := []:
      depList := []:

      if nops(basisEquations) = 0 then solutions := [seq(base[op(op(dummyList)[m])],m=1...numelems(dummyList))] else subCond := [seq(ind =~ dummyList[i],i=1..numelems(dummyList))]:

      finalEquationList := MakeUnique(map(op,[seq(`if`(subs(subCond[i],basisEquations) = [0=0],NULL,subs(subCond[i],basisEquations)),i=1...numelems(subCond))])):
      
      for i in solve(finalEquationList,maxsols=infinity) do:
      
         if lhs(i) = rhs(i) then solutions:=[op(solutions),lhs(i)] else depList :=[op(depList),lhs(i)] end if;
      
      od:
   
   end if:

   return [seq(`if`(map2(op,0,solutions[i])=base,solutions[i],NULL),i=1..nops(solutions))];
   
   end proc;

   NumberOfIndependentComponents := proc(tenss,basisEquations::list,dim::posint,startdim::integer:=0)::integer;

      description "Returns the number of independent components of a tensor, given some equations describing its behaviour and dimension";

      option remember,threadsafe;

      return nops(IndependentComponents(tenss,basisEquations,dim,startdim))

   end proc;

   Equations := proc(tenss,basisEquations::list,dim,startdim::integer:=0)::list;

   description "This takes those components that are not zero, nor are they independent, and finds expressions for them";

   option remember,threadsafe;

   uses ListTools;

      #Initialises the local variables
      local sol:=[],dummyList:=[],base,ind,rank,k,perm,i,solutions,depList,subCond,finalEquationList,deps;

      #These extract the base, indices and rank from the tensor. I can't yet find a better way to do this
      base := map2(op,0,tenss):
      ind := op(tenss):
      rank := nops([ind]):

      dummyList := [seq([seq(iquo(i,dim^r) mod dim^r,r=(rank-1)..1,-1),irem(i,dim)],i=0..(dim^(rank) - 1))];

      solutions := []:
      depList := []:

      if nops(basisEquations) = 0 then solutions := [seq(base[op(op(dummyList)[m])],m=1...numelems(dummyList))] else subCond := [seq(ind =~ dummyList[i],i=1..numelems(dummyList))]:

         finalEquationList := MakeUnique(map(op,[seq(`if`(subs(subCond[i],basisEquations) = [0=0],NULL,subs(subCond[i],basisEquations)),i=1...numelems(subCond))])):
      
         for i in solve(finalEquationList,maxsols=infinity) do:

            if lhs(i) = rhs(i) then solutions:=[op(solutions),lhs(i)] else depList :=[op(depList),lhs(i)] end if;

         od:

      end if:

      deps := op(solve(finalEquationList,depList)):
      depList := [seq(`if`(rhs(deps[i])=0,NULL,deps[i]),i=1...nops(deps))]:

      return [seq(`if`(map2(op,0,lhs(depList[i]))=base,depList[i],NULL),i=1..nops(depList))];

   end proc;

end module:
