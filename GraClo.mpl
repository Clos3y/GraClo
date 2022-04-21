GraClo := module()

   description "A module for finding the independent components and corresponding dependent components (with dependence) given a set of equations.";

   #v2.1

   # These are the functions accessible outside of the module
   export IndependentComponents,NumberOfIndependentComponents,Equations,SYM,ASYM;

   #This packages the module
   option package;

   SYM := proc(condition::list,tenss);

      #Symmetrises a tensor given a set of indices to perform the sum over

      #Initialises the local variables
      local perms := [],p,l,k,base,ind;

      #These extract the base and indices from the tensor. I can't yet find a better way to do this
      base := map2(op,0,tenss):
      ind := op(tenss):

      #Permutes the conditions
      for p in combinat:-permute(condition) do:

         #Necessary conversion for member insertion
         p := convert(p,Array):
 
         #If the index appears in the tensor, then it is ignored. Else, it is added in the correct, previous position
         for l,k in ind do:

            if member(k,condition) then else 

               p := ArrayTools:-Insert(p,l,ind[l]);
   
            end if;

         od:

      #The list of permuted indices
      perms := [op(perms),convert(p,list)];

      od:
   
   #Returns the symmetric sum
   return (1/nops(perms))*add(base[op(perms[i])],i=1...nops(perms));

   end proc:

   ASYM := proc(condition::list,tenss);

   #Antisymmetrises a tensor given a set of indices to perform the sum over
   
   #Initialises the local variables
   local perms := [],p,l,k,base,ind;

   #These extract the base and indices from the tensor. I can't yet find a better way to do this
   base := map2(op,0,tenss):
   ind := op(tenss):
   
   #Permutes the conditions
   for p in combinat:-permute(condition) do:
      
      #Necessary conversion for member insertion
      p := convert(p,Array):
      
      #If the index appears in the tensor, then it is ignored. Else, it is added in the correct, previous position
      for l,k in ind do:

         if member(k,condition) then else 

         p := ArrayTools:-Insert(p,l,ind[l]);
   
         end if;

      od:

      perms := [op(perms),convert(p,list)];

   od:

   return (1/nops(perms))*add(GroupTheory:-PermParity(Perm([seq(ListTools:-Search(s,[ind]),s in perms[i])]))*base[op(perms[i])],i=1...nops(perms));

   end proc:

   IndependentComponents := proc(tenss,basisEquations::list,dim::posint,startdim::integer:=0);

      #Find the independent components of a tensor given some equations defining its behaviour, and a dimension. One may optionally change the starting index
      
      #Initialises the local variables
      local sol:=[],dummyList:=[],base,ind,rank,k,perm,i,solutions,depList,subCond,finalEquationList;

      #These extract the base, indices and rank from the tensor. I can't yet find a better way to do this
      base := map2(op,0,tenss):
      ind := op(tenss):
      rank := nops([ind]):

      #Creates the permuted list of all numbered elements (i.e., [0,0,0],[0,0,1],[0,0,2]...)

      for k from 0 to (dim^rank - 1) do:

         perm[k] := op(ListTools:-Reverse(convert(k,'base',dim))):

         if numelems([perm[k]]) < rank then

            for i from 1 to (rank - numelems([perm[k]])) do:

               perm[k] := 0,perm[k]

            od:

         else end if:

         perm[k] := perm[k] +~ startdim;

         dummyList := [op(dummyList),[perm[k]]]:

      od: 

      #Solves the equations for the independent components

      solutions := []:
      depList := []:

      if nops(basisEquations) = 0 then solutions := [seq(base[op(op(dummyList)[m])],m=1...numelems(dummyList))] else subCond := [seq(ind =~ dummyList[i],i=1..numelems(dummyList))]:

      finalEquationList := ListTools:-MakeUnique(map(op,[seq(`if`(subs(subCond[i],basisEquations) = [0=0],NULL,subs(subCond[i],basisEquations)),i=1...numelems(subCond))])):
      
      for i in solve(finalEquationList,maxsols=infinity) do:
      
         if lhs(i) = rhs(i) then solutions:=[op(solutions),lhs(i)] else depList :=[op(depList),lhs(i)] end if;
      
      od:
   
   end if:

   return [seq(`if`(map2(op,0,solutions[i])=base,solutions[i],NULL),i=1..nops(solutions))];
   
   end proc;

   NumberOfIndependentComponents := proc(tenss,basisEquations::list,dim::posint,startdim::integer:=0)

      #Returns the number of independent components of a tensor, given some equations describing its behaviour and dimension

      return nops(IndependentComponents(tenss,basisEquations,dim,startdim))

   end proc;

   Equations := proc(tenss,basisEquations::list,dim,startdim::integer:=0);

   #This takes those components that are not zero, nor are they independent, and finds expressions for them

      #Initialises the local variables
      local sol:=[],dummyList:=[],base,ind,rank,k,perm,i,solutions,depList,subCond,finalEquationList,deps;

      #These extract the base, indices and rank from the tensor. I can't yet find a better way to do this
      base := map2(op,0,tenss):
      ind := op(tenss):
      rank := nops([ind]):

      for k from 0 to (dim^rank - 1) do:

         perm[k] := op(ListTools:-Reverse(convert(k,'base',dim))):

         if numelems([perm[k]]) < rank then

            for i from 1 to (rank - numelems([perm[k]])) do:

               perm[k] := 0,perm[k]

            od:

         else end if:

         perm[k] := perm[k] +~ startdim;

         dummyList := [op(dummyList),[perm[k]]]:
      
      od: 

      solutions := []:
      depList := []:

      if nops(basisEquations) = 0 then solutions := [seq(base[op(op(dummyList)[m])],m=1...numelems(dummyList))] else subCond := [seq(ind =~ dummyList[i],i=1..numelems(dummyList))]:

         finalEquationList := ListTools:-MakeUnique(map(op,[seq(`if`(subs(subCond[i],basisEquations) = [0=0],NULL,subs(subCond[i],basisEquations)),i=1...numelems(subCond))])):
      
         for i in solve(finalEquationList,maxsols=infinity) do:

            if lhs(i) = rhs(i) then solutions:=[op(solutions),lhs(i)] else depList :=[op(depList),lhs(i)] end if;

         od:

      end if:

      deps := op(solve(finalEquationList,depList)):
      depList := [seq(`if`(rhs(deps[i])=0,NULL,deps[i]),i=1...nops(deps))]:

      return [seq(`if`(map2(op,0,lhs(depList[i]))=base,depList[i],NULL),i=1..nops(depList))];

   end proc;

end module: