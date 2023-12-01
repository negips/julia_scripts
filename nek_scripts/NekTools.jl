#!/bin/julia

"""
  Convert from Preprocessor Coordinates indices 
  to Symmetric Coordinate indicies

  Preprocessor Corner notation  ------->  Symmetric Corner notation:

          4+-----+3    ^ s                    3+-----+4    ^ s
          /     /|     |                      /     /|     |
         /     / |     |                     /     / |     |
       8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
        |     | /     /                     |     | /     /
        |     |/     /                      |     |/     /
       5+-----+6    t                      5+-----+6    t


  NekPrepToSymm(i)

# Examples
```julia-repl
julia> j = NekPrepToSymm(3)
4
```
"""
function NekPrepToSymm(i::Int)

#  indx = [1; 2; 4; 3; 5; 6; 8; 7]

  si = i
  if i==3
    si = 4
  elseif i==4
    si = 3
  elseif i==7
    si = 8
  elseif i==8
    si = 7
  end

  if i<1 || i>8
    error("NekPrepToSymm: Index Out of bounds: $i")
  end

  return si
end  

#----------------------------------------------------------------------

"""
  Convert from Preprocessor Coordinates indices 
  to Symmetric Coordinate indicies

  Symmetric Corner notation: ------->  Preprocessor Corner notation    
                             
         3+-----+4    ^ s                     4+-----+3    ^ s                
         /     /|     |                       /     /|     |                  
        /     / |     |                      /     / |     |                  
      7+-----+8 +2    +----> r             8+-----+7 +2    +----> r           
       |     | /     /                      |     | /     /                   
       |     |/     /                       |     |/     /                    
      5+-----+6    t                       5+-----+6    t                     


  NekSymmToPrep(i)

# Examples
```julia-repl
julia> j = NekSymmToPrep(7)
8
```
"""
function NekSymmToPrep(i::Int)

#  indx = [1; 2; 4; 3; 5; 6; 8; 7]
  
  ppi = i
  if i==3
    ppi = 4
  elseif i==4
    ppi = 3
  elseif i==7
   ppi = 8
  elseif i==8
   ppi = 7
  end

  if i<1 || i>8
    error("NekSymmToPrep: Index Out of bounds: $i")
  end

  return ppi
end  
#----------------------------------------------------------------------






