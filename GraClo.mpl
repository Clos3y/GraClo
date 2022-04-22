GraClo := module()

   description "A module for finding the independent components and corresponding dependent components (with dependence) given a set of equations.":

   # These are the functions accessible outside of the module
   export IndependentComponents,NumberOfIndependentComponents,ComponentRelationships,SYM,ASYM:

   #This packages the module
   option package:

   local permuter := proc(condition::list,tensor)::Array:

      options remember,threadsafe:

      uses combinat,ArrayTools:

      #Initialises the local variables
      local indices :: exprseq := op(tensor):

      local permutations :: Array := Array(map(convert,permute(condition),Array),datatype=Array):

      #If the index appears in the tensor, then it is ignored. Else, it is added in the correct, previous position
      for local l,k in indices do:
       map[inplace](i -> `if`(member(k,condition),i,Insert(i,l,k)),permutations):
      od:

      return permutations:

   end proc:

   SYM := proc(condition::list,tensor)::`+`:

      description "Computes the symmetric sum over the specified indices":

      local i :: posint:

      options remember,threadsafe:

      uses combinat:

      local base :: symbol := map2(op,0,tensor):
     
      local permutations := permuter(condition,tensor):

      local n_perm :: posint := numbperm(condition):
         #Returns the symmetric sum
       return (1/n_perm)*add(base[op(convert(permutations[i],list))],i=1...n_perm):

   end proc:

   ASYM := proc(condition::list,tensor)::`+`:

      description "Computes the antisymmetric sum over the specified indices":

      options remember,threadsafe:

      uses combinat,GroupTheory:

      local i :: posint:

      local base :: symbol := map2(op,0,tensor):
     
      local permutations :: Array := permuter(condition,tensor):

      local n_perm :: posint := numbperm(condition):

      local signs :: list := map(PermParity,map(Perm,permute(nops(tensor)))):

     return (1/n_perm)*add(signs[i]*base[op(convert(permutations[i],list))],i=1...n_perm):

   end proc:

   local ComponentSolver := proc(tensor,dim::posint:=4,{equations := NULL,startdim::integer:=0})::list:
      option remember,threadsafe:

      local base := map2(op,0,tensor), indices := op(tensor), rank := nops(tensor):

      local i :: integer, m::posint,r::posint:

      #Creates the permuted list of all numbered elements (i.e., [0,0,0],[0,0,1],[0,0,2]...)
      local numbers_list := [seq([seq(iquo(i,dim^r) mod dim,r=(rank-1)..1,-1),irem(i,dim)],i=0..(dim^(rank) - 1))]:

      local total_perms :: posint := dim^rank:
      #Solves the equations
      local all_components :=seq(base[op(op(numbers_list)[m])],m=1...total_perms);

      if equations = NULL then:
        return [all_components]:
      else: 
      	local substituted_conditions := [seq(indices =~ numbers_list[i],i=1..total_perms)]:

	local solutions := map(x -> subs(x,equations),substituted_conditions):
	solutions := map(op,solutions):
	solutions := solve(solutions,[all_components],maxsols=infinity):
	print(select(member,R[3,2,1,0],op(solutions)));
	return op(solutions): 

      end if:

   end proc:

   IndependentComponents := proc(tensor,dim::posint:=4,{equations := NULL,startdim::integer:=0})::set:

      description "Find the independent components of a tensor given some equations defining its behaviour, and a dimension. One may optionally change the starting index":

      option remember,threadsafe:

      local solutions := ComponentSolver(tensor,dim,':-equations'=equations,':-startdim'=startdim):

      return map(lhs,select(is,solutions)):
   
   end proc:

   NumberOfIndependentComponents := proc(tensor,dim::posint:=4,{equations := NULL,startdim::integer:=0})::integer:

      description "Returns the number of independent components of a tensor, given some equations describing its behaviour and dimension":

      option remember,threadsafe:

      return nops(IndependentComponents(tensor,dim,':-equations'=equations,':-startdim'=startdim)):

   end proc:

   ComponentRelationships := proc(tensor,dim::posint:=4,{equations := NULL,startdim::integer:=0})::set:

      description "Find the independent components of a tensor given some equations defining its behaviour, and a dimension. One may optionally change the starting index":

      option remember,threadsafe:

      local solutions := ComponentSolver(tensor,dim,':-equations'=equations,':-startdim'=startdim):

      return remove(is,solutions):
   
   end proc:

end module:
